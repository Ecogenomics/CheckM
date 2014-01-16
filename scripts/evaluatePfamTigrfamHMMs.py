#!/usr/bin/env python

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

"""
Evaluate PFAM and TIGRFAM HMMs when applied to ORFs and full genomes.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os

from lib.img import IMG
from lib.pfam import PFAM
from lib.markerSet import MarkerSet

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class EvaluateHMMs(object):
    def __init__(self):
        self.pfam = PFAM()
        pass

    def translateSixFramesGenerator(self, bioseq, table=11):
        revseq = bioseq.reverse_complement()
        for i in range(3):
            yield bioseq[i:].translate(table)
            yield revseq[i:].translate(table)

    def translateSixFrames(self, genomeId):
        print '  Creating six frame translation.'
        outFile = open('./hmm_test/' + genomeId + '.six_frames.fna', 'w')

        contigFile = open(IMG.genomeDir + genomeId + '/' + genomeId + '.fna')
        for rec in SeqIO.parse(contigFile, 'fasta'):
            frameNumber = 0
            for frame in self.translateSixFramesGenerator(rec.seq):
                frameNumber += 1
                SeqIO.write(SeqRecord(frame, description='', id=rec.id + "_" + str(frameNumber)), outFile, 'fasta')

        outFile.close()

    def runPFAM(self, genomeId):
        print '  Running PFAM HMMs.'
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_ga --domtblout ./hmm_test/' + genomeId + '.pfam.table.txt ./hmm/pfam_markers.hmm ' + IMG.genomeDir + genomeId + '/' + genomeId + '.genes.derep.faa > ./hmm_test/' + genomeId + '.pfam.tsv\n')

    def runTIGRFAM(self, genomeId):
        print '  Running TIGRFAM HMMs.'
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_nc --domtblout ./hmm_test/' + genomeId + '.tigr.table.txt ./hmm/tigr_markers.hmm ' + IMG.genomeDir + genomeId + '/' + genomeId + '.genes.faa > ./hmm_test/' + genomeId + '.tigr.tsv\n')

    def runPFAM_SixFrames(self, genomeId):
        print '  Running PFAM HMMs on six frame translation.'
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_ga --domtblout ./hmm_test/' + genomeId + '.pfam.table.six_frames.txt ./hmm/pfam_markers.hmm ./hmm_test/' + genomeId + '.six_frames.fna > ./hmm_test/' + genomeId + '.pfam.six_frames.tsv\n')

    def runTIGRFAM_SixFrames(self, genomeId):
        print '  Running TIGRFAM HMMs on six frame translation.'
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_nc --domtblout ./hmm_test/' + genomeId + '.tigr.table.six_frames.txt ./hmm/tigr_markers.hmm ./hmm_test/' + genomeId + '.six_frames.fna > ./hmm_test/' + genomeId + '.tigr.six_frames.tsv\n')

    def readImgTable(self, genomeId, markers, extension, clusterIdIndex):
        hits = {}
        bHeader = True
        for line in open(IMG.genomeDir + genomeId + '/' + genomeId + extension):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            clusterId = lineSplit[clusterIdIndex]

            if clusterId in markers:
                s = hits.get(clusterId, set())
                s.add(lineSplit[0])
                hits[clusterId] = s

        return hits

        def readHmmTable(self, genomeId, extension):
            hits = {}
            for line in open('./hmm_test/' + genomeId + extension):
                if line[0] == '#' or line.strip() == '':
                    continue

                lineSplit = line.split()

                clusterId = lineSplit[4]
                if 'PF' in clusterId:
                    clusterId = clusterId.replace('PF', 'pfam')
                    clusterId = clusterId[0:clusterId.rfind('.')]

                geneId = lineSplit[0]
                evalue = float(lineSplit[12]) # i-evalue
                start = int(lineSplit[17])    # seq. start
                end = int(lineSplit[18])      # seq. end
                hits[geneId] = hits.get(geneId, []) + [[clusterId, evalue, start, end]]

            return hits

        def compareResults(self, genomeId, pfamMarkers, tigrMarkers, fout):
            # get marker hits to genes as determined by IMG
            imgPfamHits= self.readImgTable(genomeId, pfamMarkers, '.pfam.tab.txt', 8)
            imgTigrHits = self.readImgTable(genomeId, tigrMarkers, '.tigrfam.tab.txt', 6)

            print '  PFAM IMG hits: ' + str(len(imgPfamHits))
            print '  PFAM TIGR hits: ' + str(len(imgTigrHits))

            # get marker hits to genes as determined by HMMs
            hmmPfamHits = self.readHmmTable(genomeId, '.pfam.table.txt')
            hmmTigrHits = self.readHmmTable(genomeId, '.tigr.table.txt')

            print '  PFAM HMM hits: ' + str(len(hmmPfamHits))
            print '  TIGR HMM hits: ' + str(len(hmmTigrHits))

            # remove overlapping PFAM hits from the same clan
            print '  Filtering PFAM hits from the same clan.'
            filteredHmmPfamHits = self.pfam.filterHitsFromSameClan(hmmPfamHits, pfamMarkers)
            print '  Filtered PFAM hits: ' + str(len(filteredHmmPfamHits))

            # reform TIGR hits so dictionary is indexed by TIGR ids
            reformedHmmTigrHits = {}
            for geneId, hits in hmmTigrHits.iteritems():
                for h in hits:
                    tigrId = h[0]
                    s = reformedHmmTigrHits.get(tigrId, set())
                    s.add(geneId)
                    reformedHmmTigrHits[tigrId] = s

            # compare results
            pfamDiff = 0
            totalImgHits = 0
            totalHmmHits = 0
            imgAdditions = 0
            hmmAdditions = 0
            for pfamId in pfamMarkers:
                if len(imgPfamHits.get(pfamId, set())) - len (filteredHmmPfamHits.get(pfamId, set())) > 0:
                    imgAdditions += len(imgPfamHits.get(pfamId, set())) - len (filteredHmmPfamHits.get(pfamId, set()))

                if len(filteredHmmPfamHits.get(pfamId, set())) - len (imgPfamHits.get(pfamId, set())) > 0:
                    hmmAdditions += len(filteredHmmPfamHits.get(pfamId, set())) - len (imgPfamHits.get(pfamId, set()))

                pfamDiff += len(imgPfamHits.get(pfamId, set()).symmetric_difference(filteredHmmPfamHits.get(pfamId, set())))
                totalImgHits += len(imgPfamHits.get(pfamId, set()))
                totalHmmHits += len(filteredHmmPfamHits.get(pfamId, set()))

            print '  PFAM (symmetric diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(pfamDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions)
            fout.write('  PFAM (diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(pfamDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions) + '\n')

            tigrDiff = 0
            totalImgHits = 0
            totalHmmHits = 0
            imgAdditions = 0
            hmmAdditions = 0
            for tigrId in tigrMarkers:
                if len(imgTigrHits.get(tigrId, set())) - len (reformedHmmTigrHits.get(tigrId, set())) > 0:
                    imgAdditions += len(imgTigrHits.get(tigrId, set())) - len (reformedHmmTigrHits.get(tigrId, set()))

                if len(reformedHmmTigrHits.get(tigrId, set())) - len (imgTigrHits.get(tigrId, set())) > 0:
                    hmmAdditions += len(reformedHmmTigrHits.get(tigrId, set())) - len (imgTigrHits.get(tigrId, set()))

                tigrDiff += len(imgTigrHits.get(tigrId, set()).symmetric_difference(reformedHmmTigrHits.get(tigrId, set())))
                totalImgHits += len(imgTigrHits.get(tigrId, set()))
                totalHmmHits += len(reformedHmmTigrHits.get(tigrId, set()))

            print '  TIGR (symmetric diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(tigrDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions)
            print ''
            fout.write('  TIGR (diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(tigrDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions) + '\n\n')

        def compareSixFrameResults(self, genomeId, pfamMarkers, tigrMarkers, fout):
            # get marker hits to genes as determined by IMG
            imgPfamHits= self.readImgTable(genomeId, pfamMarkers, '.pfam.tab.txt', 8)
            imgTigrHits = self.readImgTable(genomeId, tigrMarkers, '.tigrfam.tab.txt', 6)

            # get marker hits to genes as determined by HMMs
            hmmPfamHits = self.readHmmTable(genomeId, pfamMarkers, '.pfam.table.six_frames.txt')
            hmmTigrHits = self.readHmmTable(genomeId, tigrMarkers, '.tigr.table.six_frames.txt')

            # compare results
            pfamDiff = 0
            totalImgHits = 0
            totalHmmHits = 0
            imgAdditions = 0
            hmmAdditions = 0
            for pfamId in pfamMarkers:
                if len(imgPfamHits.get(pfamId, set())) - len (hmmPfamHits.get(pfamId, set())) > 0:
                    imgAdditions += len(imgPfamHits.get(pfamId, set())) - len (hmmPfamHits.get(pfamId, set()))

                if len(hmmPfamHits.get(pfamId, set())) - len (imgPfamHits.get(pfamId, set())) > 0:
                    hmmAdditions += len(hmmPfamHits.get(pfamId, set())) - len (imgPfamHits.get(pfamId, set()))

                pfamDiff += abs(len(imgPfamHits.get(pfamId, set())) - len(hmmPfamHits.get(pfamId, set())))
                totalImgHits += len(imgPfamHits.get(pfamId, set()))
                totalHmmHits += len(hmmPfamHits.get(pfamId, set()))

            print '  PFAM (diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(pfamDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions)
            fout.write('  PFAM (diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(pfamDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions) + '\n')

            tigrDiff = 0
            totalImgHits = 0
            totalHmmHits = 0
            imgAdditions = 0
            hmmAdditions = 0
            for tigrId in tigrMarkers:
                if len(imgTigrHits.get(tigrId, set())) - len (hmmTigrHits.get(tigrId, set())) > 0:
                    imgAdditions += len(imgTigrHits.get(tigrId, set())) - len (hmmTigrHits.get(tigrId, set()))

                if len(hmmTigrHits.get(tigrId, set())) - len (imgTigrHits.get(tigrId, set())) > 0:
                    hmmAdditions += len(hmmTigrHits.get(tigrId, set())) - len (imgTigrHits.get(tigrId, set()))

                tigrDiff += abs(len(imgTigrHits.get(tigrId, set())) - len(hmmTigrHits.get(tigrId, set())))
                totalImgHits += len(imgTigrHits.get(tigrId, set()))
                totalHmmHits += len(hmmTigrHits.get(tigrId, set()))

            print '  TIGR (diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(tigrDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions)
            print ''
            fout.write('  TIGR (diff, IMG hits, HMM hits, IMG additional, HMM additional): ' + str(tigrDiff) + ', ' + str(totalImgHits) + ', ' + str(totalHmmHits) + ', ' + str(imgAdditions) + ', ' + str(hmmAdditions) + '\n\n')

        def run(self):
            img = IMG()

            fout = open('./data/evaluate_hmms_with_prodigal.txt', 'w', 1)

            # get list of all marker genes
            markerset = MarkerSet()
            pfamMarkers, tigrMarkers = markerset.getCalculatedMarkerGenes()

            print 'PFAM marker genes: ' + str(len(tigrMarkers))
            print 'TIGR marker genes: ' + str(len(pfamMarkers))
            print ''

            # run HMMs on each of the finished genomes
            genomeIds = img.genomeIds('Finished')
            for genomeId in genomeIds:
                print genomeId + ':'
                fout.write(genomeId + ':\n')

                self.runPFAM(genomeId)
                self.runTIGRFAM(genomeId)

                fout.write('  ORF results:\n')
                self.compareResults(genomeId, pfamMarkers, tigrMarkers, fout)

                #self.translateSixFrames(genomeId)
                #self.runPFAM_SixFrames(genomeId)
                #self.runTIGRFAM_SixFrames(genomeId)

                #fout.write('  Six-frame translation results:\n')
                #self.compareSixFrameResults(genomeId, pfamMarkers, tigrMarkers, fout)

            fout.close()

if __name__ == '__main__':
    evaluateHMMs = EvaluateHMMs()
    evaluateHMMs.run()
