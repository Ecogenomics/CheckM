###############################################################################
#
# ssuFinder.py - identify SSU (16S rRNA) in genome bins
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
import logging
import threading
import multiprocessing as mp
import concurrent.futures as cf
from collections import defaultdict

import pyhmmer

from checkm.defaultValues import DefaultValues
from checkm.common import binIdFromFilename
from checkm.util.seqUtils import readFasta, readFastaSeqIds
from checkm.util.compression import openFile


# residues in the chunks of sequences searched in parallel by nhmmer
NHMMER_CHUNK_RESIDUES = 2000000

# hits with an E-value above this default inclusion threshold of nhmmer are marked as questionable ('?')
NHMMER_INCLUSION_EVALUE = 0.01

_searchWorker = {}


class _SsuHit(object):
    """Hit of a SSU HMM to a sequence."""
    __slots__ = ('name', 'description', 'length', 'pvalue', 'evalue', 'score', 'bias', 'hmmFrom', 'hmmTo',
                 'aliFrom', 'aliTo', 'envFrom', 'envTo', 'strand')


def _initSearchWorker(hmm):
    """Set the HMM searched by a worker process."""
    _searchWorker['hmm'] = hmm
    _searchWorker['pipeline'] = pyhmmer.plan7.LongTargetsPipeline(hmm.alphabet, pyhmmer.plan7.Background(hmm.alphabet),
                                                                  E=float('inf'), incE=float('inf'))


def _searchChunk(block):
    """Search a chunk of sequences with the HMM of a worker process, returning all hits.

    The E-value of a hit is calculated by HMMER for the residues searched
    in the chunk and is returned per residue ('pvalue').
    """
    def text(value):
        return value.decode() if isinstance(value, bytes) else (value or '')

    pipeline = _searchWorker['pipeline']
    pipeline.clear()
    topHits = pipeline.search_hmm(_searchWorker['hmm'], block)

    hits = []
    for topHit in topHits:
        for d in topHit.domains:
            hit = _SsuHit()
            hit.name, hit.description, hit.length = text(topHit.name), text(topHit.description), topHit.length
            hit.pvalue = topHit.evalue / topHits.searched_residues
            hit.evalue = None
            hit.score, hit.bias, hit.strand = d.score, d.bias, d.strand
            hit.hmmFrom, hit.hmmTo = d.alignment.hmm_from, d.alignment.hmm_to
            hit.aliFrom, hit.aliTo = d.alignment.target_from, d.alignment.target_to
            hit.envFrom, hit.envTo = d.env_from, d.env_to
            hits.append(hit)

    return hits


class SSU_Finder(object):
    def __init__(self, threads):
        self.logger = logging.getLogger('timestamp')

        self.totalThreads = threads

        self.bacteriaModelFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_bacteria.hmm')
        self.archaeaModelFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_archaea.hmm')
        self.eukModelFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_euk.hmm')

        # complement of IUPAC nucleotide codes
        self.compl = str.maketrans('ACGTURYKMBDHVNacgturykmbdhvn', 'TGCAAYRMKVHDBNtgcaayrmkvhdbn')

    def _nhmmer(self, seqFile, evalue, outputPrefix, domain, modelFile):
        """Search a (possibly gzipped) sequence file with nhmmer, run through pyhmmer.

        The table of hits (nhmmer --tblout) is written to
        <outputPrefix>.<domain>.table.txt and the hits to
        <outputPrefix>.<domain>.txt, in the layout of the main output of
        nhmmer without the accuracy of each hit, which pyhmmer does not
        provide.

        pyhmmer searches the sequences with a HMM in a single thread, so
        chunks of whole sequences are searched in parallel by worker
        processes. The E-value of a hit is proportional to the number of
        residues searched, so hits are reported by workers regardless of
        their E-value, which is then calculated for all sequences. This
        gives the hits and E-values of nhmmer run with a single thread
        (nhmmer --cpu 0); with multiple threads, nhmmer splits sequences
        differently and a few hits can differ slightly.
        """
        try:
            with pyhmmer.plan7.HMMFile(modelFile) as hmmFile:
                hmms = list(hmmFile)

            with open(outputPrefix + '.' + domain + '.table.txt', 'w') as tableOut, \
                    open(outputPrefix + '.' + domain + '.txt', 'w') as fout:
                for hmm in hmms:
                    hits = self._searchChunks(hmm, seqFile, evalue)
                    self._writeHitTable(tableOut, hmm, hits, bHeader=(hmm is hmms[0]))
                    self._writeHits(fout, hmm, hits)
        except Exception as e:
            self.logger.error(f'nhmmer (pyhmmer) failed: {e!r}')
            sys.exit(1)

    def _searchChunks(self, hmm, seqFile, evalue):
        """Search chunks of sequences with a HMM in parallel, returning hits sorted by E-value.

        Sequences are read and searched in chunks of ~NHMMER_CHUNK_RESIDUES
        residues, so memory use does not depend on the size of the sequence
        file.
        """
        ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
        slots = threading.Semaphore(2 * self.totalThreads)  # limits the chunks held in memory
        totalResidues = 0
        futures = []
        with cf.ProcessPoolExecutor(self.totalThreads, mp_context=ctx, initializer=_initSearchWorker, initargs=(hmm,)) as pool, \
                openFile(seqFile, 'rb') as fin, \
                pyhmmer.easel.SequenceFile(fin, format='fasta', digital=True, alphabet=hmm.alphabet) as seqs:
            while True:
                block = seqs.read_block(residues=NHMMER_CHUNK_RESIDUES)
                if not block:
                    break

                totalResidues += sum(len(seq) for seq in block)
                slots.acquire()
                future = pool.submit(_searchChunk, block)
                future.add_done_callback(lambda _: slots.release())
                futures.append(future)

            hits = [hit for future in futures for hit in future.result()]

        # both strands are searched
        for hit in hits:
            hit.evalue = hit.pvalue * 2 * totalResidues

        hits = [hit for hit in hits if hit.evalue <= evalue]
        hits.sort(key=lambda hit: (hit.evalue, -hit.score))

        return hits

    def _writeHitTable(self, fout, hmm, hits, bHeader=True):
        """Write hits in the tabular format of 'nhmmer --tblout'."""
        queryName = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
        nameWidth = max([20] + [len(hit.name) for hit in hits])
        queryWidth = max(20, len(queryName))
        posWidth = max([7] + [len(str(hit.length)) for hit in hits])

        if bHeader:
            fout.write('#%-*s %-10s %-*s %-10s %s %s %*s %*s %*s %*s %*s %6s %9s %6s %5s  %s\n' % (
                nameWidth - 1, ' target name', 'accession', queryWidth, 'query name', 'accession', 'hmmfrom', 'hmm to',
                posWidth, 'alifrom', posWidth, 'ali to', posWidth, 'envfrom', posWidth, 'env to', posWidth, 'sq len',
                'strand', '  E-value', ' score', ' bias', 'description of target'))
            fout.write('#%*s %-10s %-*s %-10s %s %s %*s %*s %*s %*s %*s %6s %9s %6s %5s %s\n' % (
                nameWidth - 1, '-------------------', '----------', queryWidth, '--------------------', '----------',
                '-------', '-------', posWidth, '-------', posWidth, '-------', posWidth, '-------', posWidth, '-------',
                posWidth, '-------', '------', '---------', '------', '-----', '---------------------'))

        for hit in hits:
            fout.write('%-*s %-10s %-*s %-10s %7d %7d %*d %*d %*d %*d %*d %6s %9.2g %6.1f %5.1f  %s\n' % (
                nameWidth, hit.name, '-', queryWidth, queryName, '-', hit.hmmFrom, hit.hmmTo,
                posWidth, hit.aliFrom, posWidth, hit.aliTo, posWidth, hit.envFrom, posWidth, hit.envTo, posWidth, hit.length,
                '   %s  ' % hit.strand, hit.evalue, hit.score, hit.bias, hit.description or '-'))

    def _writeHits(self, fout, hmm, hits):
        """Write the hits to a HMM as in the main output of 'nhmmer --noali' (see _readHits)."""
        fout.write('Query:       %s  [M=%d]\n' % (hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name, hmm.M))
        fout.write('Annotation for each hit:\n')
        for hit in hits:
            fout.write('>> %s  %s\n' % (hit.name, hit.description))
            fout.write('   %6s %5s %9s %9s %9s %2s %9s %9s %2s %9s %9s %2s %9s\n' % (
                'score', 'bias', '  Evalue', 'hmmfrom', 'hmm to', '  ', ' alifrom ', ' ali to ', '  ',
                ' envfrom ', ' env to ', '  ', '  sq len '))
            fout.write('   %6s %5s %9s %9s %9s %2s %9s %9s %2s %9s %9s %2s %9s\n' % (
                '------', '-----', '---------', '-------', '-------', '  ', '---------', '---------', '  ',
                '---------', '---------', '  ', '---------'))
            fout.write(' %s %6.1f %5.1f %9.2g %9d %9d %s%s %9d %9d %s%s %9d %9d %s%s %9d\n\n' % (
                '!' if hit.evalue <= NHMMER_INCLUSION_EVALUE else '?', hit.score, hit.bias, hit.evalue,
                hit.hmmFrom, hit.hmmTo, '[' if hit.hmmFrom == 1 else '.', ']' if hit.hmmTo == hmm.M else '.',
                hit.aliFrom, hit.aliTo, '[' if hit.aliFrom == 1 else '.', ']' if hit.aliTo == hit.length else '.',
                hit.envFrom, hit.envTo, '[' if hit.envFrom == 1 else '.', ']' if hit.envTo == hit.length else '.',
                hit.length))
        fout.write('//\n')

    def _hmmSearch(self, seqFile, evalue, outputPrefix):
        self.logger.info('Identifying bacterial 16S.')
        self._nhmmer(seqFile, evalue, outputPrefix, 'bacteria', self.bacteriaModelFile)

        self.logger.info('Identifying archaeal 16S.')
        self._nhmmer(seqFile, evalue, outputPrefix, 'archaea', self.archaeaModelFile)

        self.logger.info('Identifying eukaryotic 18S.')
        self._nhmmer(seqFile, evalue, outputPrefix, 'euk', self.eukModelFile)

    def _readHits(self, resultsFile, domain, evalueThreshold):
        """Parse hits from nhmmer output."""
        seqInfo = {}

        bReadHit = False
        with open(resultsFile) as f:
            lines = f.readlines()

        for line in lines:
            if line[0:2] == '>>':
                lineSplit = line.split()
                seqId = lineSplit[1]
                bReadHit = True
                hitCounter = 0
                continue
            elif line.strip() == '':
                bReadHit = False
            elif bReadHit:
                hitCounter += 1
                if hitCounter >= 3:
                    lineSplit = line.split()

                    iEvalue = lineSplit[3]
                    aliFrom = int(lineSplit[7])
                    aliTo = int(lineSplit[8])

                    revComp = False
                    if aliFrom > aliTo:
                        revComp = True
                        aliFrom, aliTo = aliTo, aliFrom

                    # nhmmer coordinates are 1-based and inclusive
                    alignLen = aliTo - aliFrom + 1

                    if float(iEvalue) <= evalueThreshold:
                        seqInfo[seqId] = seqInfo.get(seqId, []) + [[domain, iEvalue, str(aliFrom), str(aliTo), str(alignLen), str(revComp)]]

        return seqInfo

    def _addHit(self, hits, seqId, info, concatenateThreshold):
        """Add hits from individual HMMs and concatenate nearby hits."""

        # check if this is the first hit to this sequence
        if seqId not in hits:
            hits[seqId] = info
            return

        # check if hits to sequence are close enough to warrant concatenating them,
        # otherwise record both hits
        baseSeqId = seqId
        index = 1
        bConcatenate = False
        concateSeqId = seqId
        while True:
            # get coordinates and strand of the new and existing hits
            startNew = int(info[2])
            endNew = int(info[3])
            bRevNew = (info[5] == 'True')

            start = int(hits[seqId][2])
            end = int(hits[seqId][3])
            bRev = (hits[seqId][5] == 'True')

            # check if hits should be concatenated
            if abs(start - endNew) < concatenateThreshold and bRevNew == bRev:
                # new hit closely preceded old hit and is on same strand
                del hits[seqId]
                info[2] = str(startNew)
                info[3] = str(end)
                info[4] = str(end - startNew + 1)
                hits[concateSeqId] = info
                bConcatenate = True

            elif abs(startNew - end) < concatenateThreshold and bRevNew == bRev:
                # new hit closely follows old hit and is on same strand
                del hits[seqId]
                info[2] = str(start)
                info[3] = str(endNew)
                info[4] = str(endNew - start + 1)
                hits[concateSeqId] = info
                bConcatenate = True

            index += 1
            newSeqId = baseSeqId + '-#' + str(index)
            if bConcatenate:
                if newSeqId in hits:
                    seqId = newSeqId  # see if other sequences concatenate
                else:
                    break
            else:
                # hits are not close enough to concatenate
                if newSeqId in hits:
                    seqId = newSeqId  # see if the new hit overlaps with this
                    concateSeqId = newSeqId
                else:
                    hits[newSeqId] = info
                    break

    def _baseSeqId(self, seqId):
        """Get sequence ID without the '-#<n>' suffix used to distinguish multiple hits to a sequence."""
        if '-#' in seqId:
            return seqId[0:seqId.rfind('-#')]

        return seqId

    def _ssuSeq(self, seq, start, end, bRevComp):
        """Extract hit from a sequence given 1-based, inclusive coordinates, reverse complementing hits on the reverse strand."""
        ssuSeq = seq[start - 1:end]
        if bRevComp:
            ssuSeq = ssuSeq.translate(self.compl)[::-1]

        return ssuSeq

    def _addDomainHit(self, hits, baseSeqId, info):
        """Add hits from different domain models and concatenate nearby hits."""
        
        startNew = int(info[2])
        endNew = int(info[3])
        lengthNew = int(info[4])

        # get all 16S hits for this contig
        curSeqIds = set()
        for curSeqId in hits:
            if self._baseSeqId(curSeqId) == baseSeqId:
                curSeqIds.add(curSeqId)

        bNewHit = True
        seqToRemove = set()
        for curSeqId in curSeqIds:
            # if hits overlap then retain only the longest
            start = int(hits[curSeqId][2])
            end = int(hits[curSeqId][3])
            length = int(hits[curSeqId][4])

            if (startNew <= start and endNew >= start) or (start <= startNew and end >= startNew):
                if lengthNew > length:
                    seqToRemove.add(curSeqId)
                else:
                    bNewHit = False
                    break

        if bNewHit:
            max_seq_num = 0
            for curSeqId in curSeqIds:
                if curSeqId in seqToRemove:
                    del hits[curSeqId]
                elif '-#' in curSeqId:
                    seq_num = int(curSeqId[curSeqId.rfind('-#') + 2:])
                    if seq_num > max_seq_num:
                        max_seq_num = seq_num

            hits['{}-#{}'.format(baseSeqId, max_seq_num+1)] = info

    def run(self, contigFile, binFiles, outputDir, evalueThreshold, concatenateThreshold):
        # make sure output directory exists
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        # get bin id of binned contigs
        self.logger.info('Determining bin assignment of sequences.')
        seqIdToBinId = {}
        for f in binFiles:
            binId = binIdFromFilename(f)
            seqIds = readFastaSeqIds(f)
            for seqId in seqIds:
                seqIdToBinId[seqId] = binId

        # identify SSU rRNA genes on contigs/scaffolds
        self.logger.info('Identifying SSU rRNAs on sequences.')
        self._hmmSearch(contigFile, evalueThreshold, os.path.join(outputDir, 'ssu'))

        # read HMM hits
        hitsPerDomain = {}
        for domain in ['archaea', 'bacteria', 'euk']:
            hits = {}

            seqInfo = self._readHits(os.path.join(outputDir, 'ssu' + '.' + domain + '.txt'), domain, evalueThreshold)
            if len(seqInfo) > 0:
                for seqId, seqHits in seqInfo.items():
                    for hit in seqHits:
                        self._addHit(hits, seqId, hit, concatenateThreshold)

            hitsPerDomain[domain] = hits

        # find best domain hit for each sequence
        bestHits = {}
        for _, hits in hitsPerDomain.items():
            for seqId, info in hits.items():
                if '-#' in seqId:
                    seqId = seqId[0:seqId.rfind('-#')]

                self._addDomainHit(bestHits, seqId, info)
                
        # relabel hits
        newBestHits = {}
        num_hits = defaultdict(int)
        for seqId, seqInfo in bestHits.items():
            if '-#' in seqId:
                seqId = seqId[0:seqId.rfind('-#')]
            
            if seqId not in num_hits:
                newBestHits[seqId] = seqInfo
            else:
                hitNum = num_hits[seqId]
                newBestHits['{}-#{}'.format(seqId, hitNum)] = seqInfo
            
            num_hits[seqId] += 1
            
        bestHits = newBestHits

        # write summary file and putative SSU rRNAs to file
        summaryFile = os.path.join(outputDir, 'ssu_summary.tsv')
        summaryOut = open(summaryFile, 'w')
        summaryOut.write('Bin Id\tSeq. Id\tHMM\ti-Evalue\tStart hit\tEnd hit\t16S/18S gene length\tRev. Complement\tSequence length\n')

        seqFile = os.path.join(outputDir, 'ssu.fna')
        seqOut = open(seqFile, 'w')

        seqs = readFasta(contigFile)

        hitsToBins = {}
        for seqId in bestHits:
            origSeqId = seqId
            if '-#' in seqId:
                seqId = seqId[0:seqId.rfind('-#')]

            if seqId in seqIdToBinId:
                binId = seqIdToBinId[seqId]
            else:
                binId = DefaultValues.UNBINNED

            seqInfo = [origSeqId] + bestHits[origSeqId]
            hitsToBins[binId] = hitsToBins.get(binId, []) + [seqInfo]

        for binId in sorted(hitsToBins.keys()):
            for seqInfo in hitsToBins[binId]:
                seqId = seqInfo[0]
                if '-#' in seqId:
                    seqId = seqId[0:seqId.rfind('-#')]

                seq = seqs[seqId]
                summaryOut.write(binId + '\t' + '\t'.join(seqInfo) + '\t' + str(len(seq)) + '\n')
                seqOut.write('>' + binId + DefaultValues.SEQ_CONCAT_CHAR + seqInfo[0] + '\n')
                seqOut.write(self._ssuSeq(seq, int(seqInfo[3]), int(seqInfo[4]), seqInfo[6] == 'True') + '\n')

        summaryOut.close()
        seqOut.close()

        self.logger.info('Identified ' + str(len(bestHits)) + ' putative SSU genes.')
        self.logger.info('Summary of identified hits written to: ' + summaryFile)
        self.logger.info('SSU sequences written to: ' + seqFile)
