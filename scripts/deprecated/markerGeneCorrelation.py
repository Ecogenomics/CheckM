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
Identify correlation between marker gene order and relative distances.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

from lib.img import IMG
from lib.markerSet import MarkerSet

from scipy.stats.mstats import spearmanr, pearsonr
from numpy import mean, std

class MarkerGeneCorelation(object):
    def __init__(self):
        pass

    def run(self):
        img = IMG()
        markerset = MarkerSet()

        print('Reading metadata.')
        metadata = img.genomeMetadata('Final')

        print('Getting marker genes.')
        pfamMarkers, tigrMarkers = markerset.getLineageMarkerGenes('Archaea')
        markerGenes = pfamMarkers.union(tigrMarkers)
        print('  Marker genes: ' + str(len(markerGenes)))

        print('Getting genomes of interest.')
        genomeIds = img.genomeIdsByTaxonomy('Archaea', 'Final')
        print('  Genomes: ' + str(len(genomeIds)))

        print('Getting position of each marker gene.')
        geneDistTable = img.geneDistTable(genomeIds, markerGenes)

        spearmanValues = []
        pearsonValues = []
        genomeIds = list(genomeIds)
        for i in range(0, len(genomeIds)):
            print(str(i+1) + ' of ' + str(len(genomeIds)))

            geneOrderI = []
            maskI = []
            for markerGenesId in markerGenes:
                if markerGenesId in geneDistTable[genomeIds[i]]:
                    geneOrderI.append(float(geneDistTable[genomeIds[i]][markerGenesId][0][0]) / metadata[genomeIds[i]]['genome size'])
                    maskI.append(0)
                else:
                    geneOrderI.append(-1)
                    maskI.append(1)


            for j in range(i+1, len(genomeIds)):
                geneOrderJ = []
                maskJ = []
                for markerGenesId in markerGenes:
                    if markerGenesId in geneDistTable[genomeIds[j]]:
                        geneOrderJ.append(float(geneDistTable[genomeIds[j]][markerGenesId][0][0]) / metadata[genomeIds[j]]['genome size'])
                        maskJ.append(0)
                    else:
                        geneOrderJ.append(-1)
                        maskJ.append(1)

                # test all translations
                bestSpearman = 0
                bestPearson = 0
                for _ in range(0, len(markerGenes)):
                    maskedI = []
                    maskedJ = []
                    for k in range(0, len(maskI)):
                        if maskI[k] == 0 and maskJ[k] == 0:
                            maskedI.append(geneOrderI[k])
                            maskedJ.append(geneOrderJ[k])
                    r, _ = spearmanr(maskedI, maskedJ)
                    if abs(r) > bestSpearman:
                        bestSpearman = abs(r)

                    r, _ = pearsonr(maskedI, maskedJ)
                    if abs(r) > bestPearson:
                        bestPearson = abs(r)

                    geneOrderJ = geneOrderJ[1:] + [geneOrderJ[0]]
                    maskJ = maskJ[1:] + [maskJ[0]]

                spearmanValues.append(bestSpearman)
                pearsonValues.append(bestPearson)

        print('Spearman: %.2f +/- %.2f: ' % (mean(spearmanValues), std(spearmanValues)))
        print('Pearson: %.2f +/- %.2f: ' % (mean(pearsonValues), std(pearsonValues)))

if __name__ == '__main__':
    markerGeneCorelation = MarkerGeneCorelation()
    markerGeneCorelation.run()
