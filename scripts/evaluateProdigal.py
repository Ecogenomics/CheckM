#!/usr/bin/env python3

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
Evaluate PFAM and TIGRFAM HMMs when applied to ORFs predicted by Prodigal and GeneMark.

[This shows that the Prodigal and GeneMark results are very similar in terms of the PFAM
  and TIGRFAM HMM hits that occur. ]
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
from checkm.pfam import PFAM
from lib.markerSet import MarkerSet

class EvaluateHMMs(object):
    def __init__(self):
        self.pfam = PFAM()
        pass

    def runProdigal(self, genomeId):
        print('  Running Prodigal.')
        os.system('time prodigal -m -c -f gff -g 11 -a ./prodigal_test/' + genomeId + '.genes.faa -d ./prodigal_test/' + genomeId + '.genes.fna -i ' + IMG.genomeDir + genomeId + '/' + genomeId + '.fna > ./prodigal_test/' + genomeId + '.prodigal.gff\n')

    def runGeneMark(self, genomeId):
        print('  Running GeneMarkS.')
        os.system('time ~/apps/genemark/gmsuite/gmsn.pl --prok --format GFF --gcode 11 --name ' + genomeId + ' --species ' + genomeId + ' --faa --fnn ' + IMG.genomeDir + genomeId + '/' + genomeId + '.fna')
        os.system('mv ' + genomeId + '.fna.faa ./genemark_test/' + genomeId + '.genes.faa')
        os.system('mv ' + genomeId + '.fna.fnn ./genemark_test/' + genomeId + '.genes.fna')
        os.system('mv ' + genomeId + '* ./genemark_test/')

    def runPFAM(self, genomeId):
        print('  Running PFAM HMMs.')
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_ga --domtblout ./genemark_test/' + genomeId + '.pfam.table.txt ./hmm/pfam_markers.hmm ./genemark_test/' + genomeId + '.genes.faa > ./genemark_test/' + genomeId + '.pfam.tsv\n')
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_ga --domtblout ./prodigal_test/' + genomeId + '.pfam.table.txt ./hmm/pfam_markers.hmm ./prodigal_test/' + genomeId + '.genes.faa > ./prodigal_test/' + genomeId + '.pfam.tsv\n')

    def runTIGRFAM(self, genomeId):
        print('  Running TIGR HMMs.')
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_nc --domtblout ./genemark_test/' + genomeId + '.tigr.table.txt ./hmm/tigr_markers.hmm ./genemark_test/' + genomeId + '.genes.faa > ./genemark_test/' + genomeId + '.tigr.tsv\n')
        os.system('hmmsearch --notextw --noali --cpu 0 --cut_nc --domtblout ./prodigal_test/' + genomeId + '.tigr.table.txt ./hmm/tigr_markers.hmm ./prodigal_test/' + genomeId + '.genes.faa > ./prodigal_test/' + genomeId + '.tigr.tsv\n')

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

        for clusterId, geneSet in list(hits.items()):
            hits[clusterId] = len(geneSet)

        return hits

    def readHmmTable(self, resultsFile):
        hits = {}
        for line in open(resultsFile):
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
        imgPfamHits = self.readImgTable(genomeId, pfamMarkers, '.pfam.tab.txt', 8)
        imgTigrHits = self.readImgTable(genomeId, tigrMarkers, '.tigrfam.tab.txt', 6)

        print(('  PFAM IMG hits: ' + str(len(imgPfamHits))))
        print(('  PFAM TIGR hits: ' + str(len(imgTigrHits))))

        # get prodigal marker hits to genes as determined by HMMs
        hmmPfamHits = self.readHmmTable('./prodigal_test/' + genomeId + '.pfam.table.txt')
        hmmTigrHits = self.readHmmTable('./prodigal_test/' + genomeId + '.tigr.table.txt')

        print(('  Prodigal PFAM HMM hits: ' + str(len(hmmPfamHits))))
        print(('  Prodigal TIGR HMM hits: ' + str(len(hmmTigrHits))))

        # remove overlapping PFAM hits from the same clan
        print('  Filtering Prodigal PFAM hits from the same clan.')
        filteredHmmPfamHits = self.pfam.filterHitsFromSameClan(hmmPfamHits, pfamMarkers)
        print(('  Filtered Prodigal PFAM hits: ' + str(len(filteredHmmPfamHits))))

        # reformat filtered PFAM hits
        prodigalPfamHits = {}
        for pfamId, geneIds in list(filteredHmmPfamHits.items()):
            prodigalPfamHits[pfamId] = len(geneIds)

        # reformat TIGR hits so dictionary is indexed by TIGR ids
        prodigalTigrHits = {}
        for _, hits in list(hmmTigrHits.items()):
            for h in hits:
                tigrId = h[0]
                prodigalTigrHits[tigrId] = prodigalTigrHits.get(tigrId, 0) + 1

        # get GeneMark marker hits to genes as determined by HMMs
        hmmPfamHits = self.readHmmTable('./genemark_test/' + genomeId + '.pfam.table.txt')
        hmmTigrHits = self.readHmmTable('./genemark_test/' + genomeId + '.tigr.table.txt')

        print(('  GeneMark PFAM HMM hits: ' + str(len(hmmPfamHits))))
        print(('  GeneMark TIGR HMM hits: ' + str(len(hmmTigrHits))))

        # remove overlapping PFAM hits from the same clan
        print('  Filtering GeneMark PFAM hits from the same clan.')
        filteredHmmPfamHits = self.pfam.filterHitsFromSameClan(hmmPfamHits, pfamMarkers)
        print(('  Filtered GeneMark PFAM hits: ' + str(len(filteredHmmPfamHits))))

        # reformat filtered PFAM hits
        genemarkPfamHits = {}
        for pfamId, geneIds in list(filteredHmmPfamHits.items()):
            genemarkPfamHits[pfamId] = len(geneIds)

        # reformat TIGR hits so dictionary is indexed by TIGR ids
        genemarkTigrHits = {}
        for _, hits in list(hmmTigrHits.items()):
            for h in hits:
                tigrId = h[0]
                genemarkTigrHits[tigrId] = genemarkTigrHits.get(tigrId, 0) + 1

        # compare results
        prodigalDiff = 0
        genemarkDiff = 0
        totalImgHits = 0
        totalProdigalHits = 0
        totalGeneMarkHits = 0
        for pfamId in pfamMarkers:
            prodigalDiff += abs(prodigalPfamHits.get(pfamId, 0) - imgPfamHits.get(pfamId, 0))
            genemarkDiff += abs(genemarkPfamHits.get(pfamId, 0) - imgPfamHits.get(pfamId, 0))
            totalImgHits += imgPfamHits.get(pfamId, 0)
            totalProdigalHits += prodigalPfamHits.get(pfamId, 0)
            totalGeneMarkHits += genemarkPfamHits.get(pfamId, 0)

        print(('  PFAM (Prodigal diff, GeneMark diff, IMG hits, Prodigal hits, GeneMark hits): ' + str(prodigalDiff) + ', ' + str(genemarkDiff) + ', ' + str(totalImgHits) + ', ' + str(totalProdigalHits) + ', ' + str(totalGeneMarkHits)))
        fout.write('  PFAM (Prodigal diff, GeneMark diff, IMG hits, Prodigal hits, GeneMark hits): ' + str(prodigalDiff) + ', ' + str(genemarkDiff) + ', ' + str(totalImgHits) + ', ' + str(totalProdigalHits) + ', ' + str(totalGeneMarkHits) + '\n')

        prodigalDiff = 0
        genemarkDiff = 0
        totalImgHits = 0
        totalProdigalHits = 0
        totalGeneMarkHits = 0
        for tigrId in tigrMarkers:
            prodigalDiff += abs(prodigalTigrHits.get(tigrId, 0) - imgTigrHits.get(tigrId, 0))
            genemarkDiff += abs(genemarkTigrHits.get(tigrId, 0) - imgTigrHits.get(tigrId, 0))
            totalImgHits += imgTigrHits.get(tigrId, 0)
            totalProdigalHits += prodigalTigrHits.get(tigrId, 0)
            totalGeneMarkHits += genemarkTigrHits.get(tigrId, 0)

        print(('  TIGR (Prodigal diff, GeneMark diff, IMG hits, Prodigal hits, GeneMark hits): ' + str(prodigalDiff) + ', ' + str(genemarkDiff) + ', ' + str(totalImgHits) + ', ' + str(totalProdigalHits) + ', ' + str(totalGeneMarkHits)))
        print('')
        fout.write('  TIGR (Prodigal diff, GeneMark diff, IMG hits, Prodigal hits, GeneMark hits): ' + str(prodigalDiff) + ', ' + str(genemarkDiff) + ', ' + str(totalImgHits) + ', ' + str(totalProdigalHits) + ', ' + str(totalGeneMarkHits) + '\n\n')

    def run(self):
        img = IMG()

        fout = open('./data/evaluate_prodigal.txt', 'w', 1)

        # get list of all marker genes
        markerset = MarkerSet()
        pfamMarkers, tigrMarkers = markerset.getCalculatedMarkerGenes()

        print(('PFAM marker genes: ' + str(len(tigrMarkers))))
        print(('TIGR marker genes: ' + str(len(pfamMarkers))))
        print('')

        # run HMMs on each of the finished genomes
        genomeIds = img.genomeIds('Finished')
        for genomeId in genomeIds:
            print((genomeId + ':'))
            fout.write(genomeId + ':\n')

            self.runProdigal(genomeId)
            self.runGeneMark(genomeId)

            self.runPFAM(genomeId)
            self.runTIGRFAM(genomeId)

            self.compareResults(genomeId, pfamMarkers, tigrMarkers, fout)

        fout.close()

if __name__ == '__main__':
    evaluateHMMs = EvaluateHMMs()
    evaluateHMMs.run()
