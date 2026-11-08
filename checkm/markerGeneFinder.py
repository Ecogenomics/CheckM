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

import gzip
import os
import sys
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
from checkm.util.pfam import PFAM


def _fetchHmms(accessions, hmmModelFile):
    """Create a HMM file with the specified models and return the parsed models.

    Run in a worker process.
    """
    MarkerSetParser().fetchHmms(accessions, hmmModelFile)
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
        aaGeneFile = prodigal.aaGeneFile
    else:
        aaGeneFile = binFile
        aaGeneFileSave = os.path.join(binDir, DefaultValues.PRODIGAL_AA)
        if not aaGeneFile.endswith(".gz"):
            shutil.copyfile(aaGeneFile, aaGeneFileSave)
        else:
            with gzip.open(aaGeneFile, 'rt') as f_in:
                with open(aaGeneFileSave, 'w') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        aaGeneFile = aaGeneFileSave

    # skip bins without any genes, as there is nothing to search
    if os.stat(aaGeneFile).st_size == 0:
        return False

    # run HMMER
    keepAlignStr = ''
    if not bKeepAlignment:
        keepAlignStr = '--noali'

    hmmer = HMMERRunner()
    hmmer.search(hmmModelFile, aaGeneFile,
                 os.path.join(binDir, tableOut), os.path.join(binDir, hmmerOut),
                 '--cpu ' + str(threadsPerSearch) + ' --notextw -E 0.1 --domE 0.1 ' + keepAlignStr,
                 bKeepAlignment)

    return True


def _fileSize(f):
    """Size of a file in bytes, or 0 if it can't be determined."""
    try:
        return os.path.getsize(f)
    except OSError:
        return 0


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

        # make sure HMMER and prodigal are on system path
        HMMERRunner()

        if not bCalledGenes:
            ProdigalRunner('')

        self.threadsPerSearch = max(1, int(self.totalThreads / len(binFiles)))
        self.logger.info("Identifying marker genes in %d bins with %d threads:" % (len(binFiles), self.totalThreads))

        # determine the HMMs to search for in each bin, grouping bins searched with the same HMMs
        groups = self._groupBinsByHmms(binFiles, markerFile)

        tmpDir = tempfile.mkdtemp(prefix='checkm_hmms_')
        try:
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

        pfam = PFAM(DefaultValues.PFAM_CLAN_FILE)
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
        alone at the end.
        """

        # order bins and groups from largest to smallest
        orderedGroups = []
        for accessions, groupBinFiles in groups:
            groupBinFiles = sorted(groupBinFiles, key=lambda f: (-_fileSize(f), f))
            orderedGroups.append((accessions, groupBinFiles))
        orderedGroups.sort(key=lambda g: (-sum(_fileSize(f) for f in g[1]), g[1][0]))

        # models of a HMM marker file are parsed once
        markerFileModels = None
        if any(accessions is None for accessions, _ in orderedGroups):
            markerFileModels = HmmModelParser(markerFile).models()

        numBins = sum(len(groupBinFiles) for _, groupBinFiles in orderedGroups)
        numProcessedBins = 0
        self._reportProgress(numProcessedBins, numBins)

        binIdToModels = {}
        skippedBinIds = []

        pendingGroups = deque(range(len(orderedGroups)))  # groups not yet started
        readySearches = deque()      # (group index, bin file) of bins whose HMM file is ready
        groupHmmFile = {}            # HMM file of each started group
        groupModels = {}             # HMM models of each group with a ready HMM file
        groupRemaining = {}          # number of bins in a group not yet searched
        running = {}                 # future -> (task, group index, bin file)

        ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
        pool = cf.ProcessPoolExecutor(max_workers=self.totalThreads, mp_context=ctx)
        try:
            while pendingGroups or readySearches or running:
                # keep all workers busy, preferring searches of bins whose HMM file is ready
                while len(running) < self.totalThreads:
                    if readySearches:
                        groupIdx, binFile = readySearches.popleft()
                        future = pool.submit(_searchBin, binFile, groupHmmFile[groupIdx], outDir,
                                             tableOut, hmmerOut, bKeepAlignment, bNucORFs,
                                             bCalledGenes, self.threadsPerSearch)
                        running[future] = ('search', groupIdx, binFile)
                    elif pendingGroups and len(groupHmmFile) < self.totalThreads:
                        groupIdx = pendingGroups.popleft()
                        accessions, groupBinFiles = orderedGroups[groupIdx]
                        groupRemaining[groupIdx] = len(groupBinFiles)
                        if accessions is None:
                            # search with the marker file directly
                            groupHmmFile[groupIdx] = markerFile
                            groupModels[groupIdx] = markerFileModels
                            readySearches.extend((groupIdx, f) for f in groupBinFiles)
                        else:
                            hmmModelFile = os.path.join(tmpDir, f'markers_{groupIdx}.hmm')
                            groupHmmFile[groupIdx] = hmmModelFile
                            future = pool.submit(_fetchHmms, accessions, hmmModelFile)
                            running[future] = ('fetch', groupIdx, None)
                    else:
                        break

                done, _ = cf.wait(list(running), return_when=cf.FIRST_COMPLETED)
                for future in done:
                    task, groupIdx, binFile = running.pop(future)
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
                        groupModels[groupIdx] = result
                        readySearches.extend((groupIdx, f) for f in orderedGroups[groupIdx][1])
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

                    numProcessedBins += 1
                    self._reportProgress(numProcessedBins, numBins)
        except BaseException:
            # stop remaining work (e.g., after a failure or keyboard interrupt)
            self._terminatePool(pool)
            raise

        pool.shutdown()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

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

    def _reportProgress(self, numProcessedBins, numBins):
        """Report number of processed bins."""
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()
