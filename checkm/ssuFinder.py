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
import logging
from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.common import binIdFromFilename
from checkm.util.seqUtils import readFasta, readFastaSeqIds


class SSU_Finder(object):
    def __init__(self, threads):
        self.logger = logging.getLogger('timestamp')

        self.totalThreads = threads

        self.bacteriaModelFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_bacteria.hmm')
        self.archaeaModelFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_archaea.hmm')
        self.eukModelFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_euk.hmm')

    def _hmmSearch(self, seqFile, evalue, outputPrefix):
        if seqFile.endswith('gz'):
            pipe = 'zcat ' + seqFile + ' | '
        else:
            pipe = 'cat ' + seqFile + ' | '

        self.logger.info('Identifying bacterial 16S.')
        os.system(pipe + 'nhmmer --noali --cpu ' + str(self.totalThreads) + ' -o ' + outputPrefix + '.bacteria.txt --tblout ' + outputPrefix + '.bacteria.table.txt -E ' + str(evalue) + ' ' + self.bacteriaModelFile + ' -')

        self.logger.info('Identifying archaeal 16S.')
        os.system(pipe + 'nhmmer --noali --cpu ' + str(self.totalThreads) + ' -o ' + outputPrefix + '.archaea.txt --tblout ' + outputPrefix + '.archaea.table.txt -E ' + str(evalue) + ' ' + self.archaeaModelFile + ' -')

        self.logger.info('Identifying eukaryotic 18S.')
        os.system(pipe + 'nhmmer --noali --cpu ' + str(self.totalThreads) + ' -o ' + outputPrefix + '.euk.txt --tblout ' + outputPrefix + '.euk.table.txt -E ' + str(evalue) + ' ' + self.eukModelFile + ' -')

    def _readHits(self, resultsFile, domain, evalueThreshold):
        """Parse hits from nhmmer output."""
        seqInfo = {}

        bReadHit = False
        for line in open(resultsFile):
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

                    alignLen = int(aliTo) - int(aliFrom)

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
            # if hits overlap then retain only the longest
            startNew = int(info[2])
            endNew = int(info[3])
            bRevNew = bool(info[5])

            start = int(hits[seqId][2])
            end = int(hits[seqId][3])
            bRev = bool(info[5])

            # check if hits should be concatenated
            if abs(start - endNew) < concatenateThreshold and bRevNew == bRev:
                # new hit closely preceded old hit and is on same strand
                del hits[seqId]
                info[2] = str(startNew)
                info[3] = str(end)
                info[4] = str(end - startNew)
                hits[concateSeqId] = info
                bConcatenate = True

            elif abs(startNew - end) < concatenateThreshold and bRevNew == bRev:
                # new hit closely follows old hit and is on same strand
                del hits[seqId]
                info[2] = str(start)
                info[3] = str(endNew)
                info[4] = str(endNew - start)
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

    def _addDomainHit(self, hits, baseSeqId, info):
        """Add hits from different domain models and concatenate nearby hits."""
        
        startNew = int(info[2])
        endNew = int(info[3])
        lengthNew = int(info[4])

        # get all 16S hits for this contig
        curSeqIds = set()
        for curSeqId in hits:
            if baseSeqId in curSeqId:
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
                    seq_num = int(curSeqId.split('-#')[1])
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

        # identify 16S reads from contigs/scaffolds
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
                seqOut.write(seq[int(seqInfo[3]) - 1:int(seqInfo[4]) - 1] + '\n')

        summaryOut.close()
        seqOut.close()

        self.logger.info('Identified ' + str(len(bestHits)) + ' putative SSU genes.')
        self.logger.info('Summary of identified hits written to: ' + summaryFile)
        self.logger.info('SSU sequences written to: ' + seqFile)
