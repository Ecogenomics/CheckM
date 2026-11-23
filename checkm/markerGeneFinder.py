###############################################################################
#
# markerGeneFinder.py - identify marker genes in genome bins
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import bisect
import shutil
import logging
import tempfile
import multiprocessing as mp
import concurrent.futures as cf
from collections import deque

from checkm.common import binIdFromFilename, makeSurePathExists
from checkm.defaultValues import DefaultValues

from checkm.hmmer import HMMERRunner
from checkm.prodigal import ProdigalRunner

from checkm.markerSets import MarkerSetParser, BinMarkerSets
from checkm.hmmerModelParser import HmmModelParser
from checkm.util.pfam import sharedPfam
from checkm.util.forking import frozenHeap, runInChild
from checkm.util.progress import progressBar
from checkm.util.compression import compressFile, hasContent, resolvePath


def _fetchHmms(accessions, hmmModelFile):
    """Create a HMM file with the specified models and return the parsed models.

    Run in a worker process.
    """
    runInChild(MarkerSetParser().fetchHmms, accessions, hmmModelFile)
    return HmmModelParser(hmmModelFile).models()


def _searchBin(binFile, hmmModelFile, outDir, tableOut, hmmerOut, bKeepAlignment, bNucORFs, bCalledGenes, threadsPerSearch):
    """Call genes (unless provided) and search them with HMMER.

    Run in a worker process. Returns False if the bin contains no genes,
    in which case no search is performed, otherwise True.
    """

    binId = binIdFromFilename(binFile)
    binDir = os.path.join(outDir, 'bins', binId)
    makeSurePathExists(binDir)

    # run Prodigal
    if not bCalledGenes:
        prodigal = ProdigalRunner(binDir)
        if not prodigal.areORFsCalled(bNucORFs):
            prodigal.run(binFile, bNucORFs)
        aaGeneFile = resolvePath(prodigal.aaGeneFile)  # genes may be from CheckM <= 1.2.x (genes.faa)
    else:
        # save a gzip-compressed copy of the genes (which may be gzip-compressed)
        aaGeneFile = compressFile(binFile, os.path.join(binDir, DefaultValues.PRODIGAL_AA), remove=False)

    # skip bins without any genes, as there is nothing to search
    if not hasContent(aaGeneFile):
        return False

    # run HMMER (through pyhmmer)
    runInChild(HMMERRunner().search, hmmModelFile, aaGeneFile,
                os.path.join(binDir, tableOut), os.path.join(binDir, hmmerOut),
                evalue=0.1, domEvalue=0.1, bKeepOutput=bKeepAlignment, threads=threadsPerSearch)

    return True


def _fileSize(f):
    """Size of a file in bytes, or 0 if it can't be determined."""
    try:
        return os.path.getsize(f)
    except OSError:
        return 0


def _genomeSize(f):
    """Estimated size of a genome file in bytes.

    This is the file size or, for a gzip-compressed file, the uncompressed
    size recorded in its trailer (modulo 4 GiB).
    """
    size = _fileSize(f)
    if f.endswith('.gz') and size >= 4:
        try:
            with open(f, 'rb') as fh:
                fh.seek(-4, os.SEEK_END)
                return int.from_bytes(fh.read(4), 'little')
        except OSError:
            pass

    return size


# Prodigal's peak memory is ~28 bytes per base of a genome plus ~10 MB (Prodigal 2.6.3;
# genomes of 0.3 to 10 Mbp), and its run time is roughly proportional to genome size
PRODIGAL_MEMORY_PER_BASE = 28
PRODIGAL_MEMORY_BASE = 10 * 1024**2

# memory budget for running Prodigal in parallel, relative to the mean memory of the shortest possible run
GENE_CALLING_MEMORY_FACTOR = 1.2


def _prodigalMemory(genomeSize):
    """Estimated peak memory of Prodigal for a genome, in bytes."""
    return PRODIGAL_MEMORY_PER_BASE * genomeSize + PRODIGAL_MEMORY_BASE


def _geneCallingMemoryBudget(genomeSizes, threads):
    """Memory budget for calling genes in the given genomes in parallel, in bytes.

    Both the memory and the run time of Prodigal grow with genome size, so
    the memory-time product of a set of genomes is fixed, whatever the order
    in which they are processed. The budget is a little more than the mean
    memory of the shortest possible run (the total work divided by the
    number of threads, or the time of the largest genome), and at least the
    memory of the largest genome. Starting the largest genomes first without
    a budget instead runs the largest genomes all at once.
    """
    if not genomeSizes:
        return 0

    makespan = max(sum(genomeSizes) / threads, max(genomeSizes))
    meanMemory = sum(_prodigalMemory(s) * s for s in genomeSizes) / makespan

    return max(GENE_CALLING_MEMORY_FACTOR * meanMemory, _prodigalMemory(max(genomeSizes)))


def _needsGeneCalling(binFile, outDir, bNucORFs, bCalledGenes):
    """Check if genes must be called for a bin (see _searchBin and ProdigalRunner.areORFsCalled)."""
    if bCalledGenes:
        return False

    binDir = os.path.join(outDir, 'bins', binIdFromFilename(binFile))
    geneFile = DefaultValues.PRODIGAL_NT if bNucORFs else DefaultValues.PRODIGAL_AA

    return not hasContent(os.path.join(binDir, geneFile))


class _ReadyBins():
    """Bins ready to be searched, started so that the memory of Prodigal stays within a budget.

    Bins are started in the order in which their groups became ready and,
    within a group, from largest to smallest, so a few large bins are not
    left running alone at the end. A bin whose gene calling does not fit
    within the remaining budget is passed over for the largest bin that
    does; if no bin fits, none is started until memory is released. As the
    budget is at least the memory of the largest bin, a bin can always be
    started when no genes are being called. Bins whose genes have already
    been called need no memory for gene calling and always fit.
    """

    def __init__(self, budget):
        self.budget = budget
        self.memoryInUse = 0
        self.groups = {}  # group index -> [(estimated memory, -rank, bin file)] in ascending order

    def __bool__(self):
        return bool(self.groups)

    def add(self, groupIdx, binFiles, geneCallingMemory):
        """Add the bins of a group, given from largest to smallest."""
        if binFiles:
            self.groups[groupIdx] = sorted((geneCallingMemory[f], -rank, f) for rank, f in enumerate(binFiles))

    def pop(self):
        """Start the next bin, returning (group index, bin file, memory), or None if no bin fits."""
        room = self.budget - self.memoryInUse
        for groupIdx, bins in self.groups.items():
            # the largest bin that fits (the first in rank among bins of equal memory)
            idx = bisect.bisect_right(bins, room, key=lambda b: b[0]) - 1
            if idx >= 0:
                memory, _, binFile = bins.pop(idx)
                if not bins:
                    del self.groups[groupIdx]
                self.memoryInUse += memory
                return groupIdx, binFile, memory

        return None

    def release(self, memory):
        """Release the memory of a bin once it has been searched."""
        self.memoryInUse -= memory


class MarkerGeneFinder():
    """Identify marker genes within binned sequences using Prodigal and HMMER."""
    def __init__(self, threads):
        self.logger = logging.getLogger('timestamp')
        self.totalThreads = threads

    def find(self, binFiles, outDir, tableOut, hmmerOut, markerFile, bKeepAlignment, bNucORFs, bCalledGenes):
        """Identify marker genes in each bin using prodigal and HMMER.

        Returns a dictionary mapping bin IDs to HMM models, sorted by
        bin ID. Bins without any genes are skipped with a warning and
        are not included in the returned dictionary. Bins that are
        searched with the same HMMs share the same model dictionary.
        """

        # make sure prodigal is on system path
        if not bCalledGenes:
            ProdigalRunner('')

        self.threadsPerSearch = max(1, int(self.totalThreads / len(binFiles)))
        self.logger.info("Identifying marker genes in %d bins with %d threads:" % (len(binFiles), self.totalThreads))

        # determine the HMMs to search for in each bin, grouping bins searched with the same HMMs
        groups = self._groupBinsByHmms(binFiles, markerFile)

        tmpDir = tempfile.mkdtemp(prefix='checkm_hmms_')
        try:
            # the heap is frozen while worker processes may be started (the pool starts them when
            # work is first submitted), so a garbage collection in a worker does not copy the memory
            # of this process (e.g., ~0.3 GiB of marker sets in analyze for 1,000 bins, which
            # otherwise grew to ~10 GiB with 64 workers)
            with frozenHeap():
                binIdToModels, skippedBinIds = self._processGroups(groups, markerFile, tmpDir, outDir,
                                                                   tableOut, hmmerOut, bKeepAlignment,
                                                                   bNucORFs, bCalledGenes)
        finally:
            shutil.rmtree(tmpDir, ignore_errors=True)

        for binId in sorted(skippedBinIds):
            if bCalledGenes:
                self.logger.warning(f'Skipping bin {binId}: the gene file contains no genes.')
            else:
                self.logger.warning(f'Skipping bin {binId}: no genes were identified by Prodigal.')

            # remove the bin's (empty) output directory so commands that identify
            # bins from the output directory (e.g., tree_qa, lineage_set, qa) skip it too
            shutil.rmtree(os.path.join(outDir, 'bins', binId), ignore_errors=True)

        # bins in a deterministic order
        return {binId: binIdToModels[binId] for binId in sorted(binIdToModels)}

    def _groupBinsByHmms(self, binFiles, markerFile):
        """Group bins by the HMMs they must be searched with.

        Returns a list of (HMM accessions, bin files) tuples. The HMM
        accessions are None if bins are searched with all HMMs in the
        marker file (i.e., the marker file is a HMM file). The marker
        file is parsed once, rather than once per bin.
        """

        markerSetParser = MarkerSetParser()
        markerFileType = markerSetParser.markerFileType(markerFile)
        if markerFileType == BinMarkerSets.HMM_MODELS_SET:
            return [(None, list(binFiles))]

        if markerFileType == BinMarkerSets.TAXONOMIC_MARKER_SET:
            binMarkerSets = markerSetParser.parseTaxonomicMarkerSetFile(markerFile)
            binIdToBinMarkerSets = {binIdFromFilename(f): binMarkerSets for f in binFiles}
        else:
            binIdToBinMarkerSets = markerSetParser.parseLineageMarkerSetFile(markerFile)

        markerGenesToBinFiles = {}
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)
            if binId not in binIdToBinMarkerSets:
                self.logger.error(f'Bin {binId} is not in the marker file {markerFile}.')
                sys.exit(1)

            markerGenes = frozenset(binIdToBinMarkerSets[binId].getMarkerGenes())
            markerGenesToBinFiles.setdefault(markerGenes, []).append(binFile)

        pfam = sharedPfam(DefaultValues.PFAM_CLAN_FILE)
        groups = []
        for markerGenes, groupBinFiles in markerGenesToBinFiles.items():
            accessions = frozenset(markerGenes | pfam.genesInSameClan(markerGenes))
            groups.append((accessions, groupBinFiles))

        return groups

    def _processGroups(self, groups, markerFile, tmpDir, outDir, tableOut, hmmerOut, bKeepAlignment, bNucORFs, bCalledGenes):
        """Search bins in parallel, creating the HMM file of each group once.

        A group's HMM file is created by a worker when first needed and
        removed once all bins in the group have been searched. At most one
        HMM file per thread exists at a time to bound disk usage. Larger
        bins are processed first so a few large bins are not left running
        alone at the end, but bins whose genes must be called are started
        only if the estimated memory of Prodigal fits within a budget (see
        _ReadyBins), so the largest genomes are not all called at once.
        """

        # order bins and groups from largest to smallest
        genomeSize = {f: _genomeSize(f) for _, groupBinFiles in groups for f in groupBinFiles}
        orderedGroups = []
        for accessions, groupBinFiles in groups:
            groupBinFiles = sorted(groupBinFiles, key=lambda f: (-genomeSize[f], f))
            orderedGroups.append((accessions, groupBinFiles))
        orderedGroups.sort(key=lambda g: (-sum(genomeSize[f] for f in g[1]), g[1][0]))

        # estimated memory of Prodigal for bins whose genes must be called (0 for other bins)
        geneCallingMemory = {}
        geneCallingSizes = []
        for _, groupBinFiles in orderedGroups:
            for f in groupBinFiles:
                if _needsGeneCalling(f, outDir, bNucORFs, bCalledGenes):
                    geneCallingMemory[f] = _prodigalMemory(genomeSize[f])
                    geneCallingSizes.append(genomeSize[f])
                else:
                    geneCallingMemory[f] = 0

        # models of a HMM marker file are parsed once
        markerFileModels = None
        if any(accessions is None for accessions, _ in orderedGroups):
            markerFileModels = HmmModelParser(markerFile).models()

        numBins = sum(len(groupBinFiles) for _, groupBinFiles in orderedGroups)
        progress = progressBar(numBins, 'Processing bins', 'bin')

        binIdToModels = {}
        skippedBinIds = []

        pendingGroups = deque(range(len(orderedGroups)))  # groups not yet started
        readySearches = _ReadyBins(_geneCallingMemoryBudget(geneCallingSizes, self.totalThreads))  # bins whose HMM file is ready
        groupHmmFile = {}            # HMM file of each started group
        groupModels = {}             # HMM models of each group with a ready HMM file
        sharedModels = {}            # HMM accession -> model; groups share the models of HMMs they have in common
        groupRemaining = {}          # number of bins in a group not yet searched
        running = {}                 # future -> (task, group index, bin file, estimated memory of gene calling)

        ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
        pool = cf.ProcessPoolExecutor(max_workers=self.totalThreads, mp_context=ctx)
        try:
            while pendingGroups or readySearches or running:
                # keep all workers busy, preferring searches of bins whose HMM file is ready
                while len(running) < self.totalThreads:
                    nextSearch = readySearches.pop() if readySearches else None
                    if nextSearch:
                        groupIdx, binFile, memory = nextSearch
                        future = pool.submit(_searchBin, binFile, groupHmmFile[groupIdx], outDir,
                                             tableOut, hmmerOut, bKeepAlignment, bNucORFs,
                                             bCalledGenes, self.threadsPerSearch)
                        running[future] = ('search', groupIdx, binFile, memory)
                    elif pendingGroups and len(groupHmmFile) < self.totalThreads:
                        groupIdx = pendingGroups.popleft()
                        accessions, groupBinFiles = orderedGroups[groupIdx]
                        groupRemaining[groupIdx] = len(groupBinFiles)
                        if accessions is None:
                            # search with the marker file directly
                            groupHmmFile[groupIdx] = markerFile
                            groupModels[groupIdx] = markerFileModels
                            readySearches.add(groupIdx, groupBinFiles, geneCallingMemory)
                        else:
                            hmmModelFile = os.path.join(tmpDir, f'markers_{groupIdx}.hmm')
                            groupHmmFile[groupIdx] = hmmModelFile
                            future = pool.submit(_fetchHmms, accessions, hmmModelFile)
                            running[future] = ('fetch', groupIdx, None, 0)
                    else:
                        break

                done, _ = cf.wait(list(running), return_when=cf.FIRST_COMPLETED)
                for future in done:
                    task, groupIdx, binFile, memory = running.pop(future)
                    readySearches.release(memory)
                    try:
                        result = future.result()  # re-raises any error from the worker
                    except (Exception, SystemExit) as e:
                        failedItem = f'bin {binIdFromFilename(binFile)}' if task == 'search' else 'the HMM file of a marker set'
                        if isinstance(e, (SystemExit, cf.process.BrokenProcessPool)):
                            # the worker has already reported the error (or was killed)
                            self.logger.error(f'Failed to process {failedItem}. See the error messages above.')
                        else:
                            self.logger.error(f'Failed to process {failedItem}: {e!r}')
                        sys.exit(1)

                    if task == 'fetch':
                        # all models are from the same marker file, so a model parsed for an earlier
                        # group is reused; memory then depends on the number of HMMs, not groups
                        groupModels[groupIdx] = {acc: sharedModels.setdefault(acc, model) for acc, model in result.items()}
                        readySearches.add(groupIdx, orderedGroups[groupIdx][1], geneCallingMemory)
                        continue

                    binId = binIdFromFilename(binFile)
                    if result:
                        binIdToModels[binId] = groupModels[groupIdx]
                    else:
                        skippedBinIds.append(binId)

                    groupRemaining[groupIdx] -= 1
                    if groupRemaining[groupIdx] == 0:
                        # all bins in group searched
                        if orderedGroups[groupIdx][0] is not None:
                            os.remove(groupHmmFile[groupIdx])
                        del groupHmmFile[groupIdx]
                        del groupModels[groupIdx]

                    progress.update(1)
        except BaseException:
            # stop remaining work (e.g., after a failure or keyboard interrupt)
            self._terminatePool(pool)
            raise
        finally:
            progress.close()

        pool.shutdown()

        return binIdToModels, skippedBinIds

    def _terminatePool(self, pool):
        """Stop all worker processes of a process pool."""
        processes = getattr(pool, '_processes', None) or {}
        for p in list(processes.values()):
            try:
                p.terminate()
            except Exception:
                pass
        pool.shutdown(wait=False)
