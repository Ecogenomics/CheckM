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
Calculate size of marker set.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse
import random

from lib.markerSet import MarkerSet
from lib.img import IMG
from lib.plots.heatmap import Heatmap

class PlotMarkerSetDistribution(object):
    def __init__(self):
        pass

    def run(self, taxonomyStr, ubiquityThreshold, singleCopyThreshold, numBins, numRndGenomes):
        img = IMG()
        markerSet = MarkerSet()

        metadata = img.genomeMetadata()
        lineageGenomeIds = img.genomeIdsByTaxonomy(taxonomyStr, metadata)

        # build marker set from finished prokaryotic genomes
        genomeIds = []
        for genomeId in lineageGenomeIds:
            if metadata[genomeId]['status'] == 'Finished' and (metadata[genomeId]['taxonomy'][0] == 'Bacteria' or metadata[genomeId]['taxonomy'][0] == 'Archaea'):
                genomeIds.append(genomeId)
        genomeIds = set(genomeIds) - img.genomesWithMissingData(genomeIds)

        print(('Lineage ' + taxonomyStr + ' contains ' + str(len(genomeIds)) + ' genomes.'))

        # get marker set
        countTable = img.countTable(genomeIds)
        countTable = img.filterTable(genomeIds, countTable, 0.9*ubiquityThreshold, 0.9*singleCopyThreshold)
        markerGenes = markerSet.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
        tigrToRemove = img.identifyRedundantTIGRFAMs(markerGenes)
        markerGenes = markerGenes - tigrToRemove
        geneDistTable = img.geneDistTable(genomeIds, markerGenes, spacingBetweenContigs=1e6)

        print(('Number of marker genes: ' + str(len(markerGenes))))

        # randomly set genomes to plot
        if numRndGenomes != -1:
            genomeIds = random.sample(list(genomeIds), numRndGenomes)
        genomeIds = set(genomeIds)

        # plot distribution of marker genes
        filename = 'geneDistribution.' + taxonomyStr.replace(';','_') + '.' + str(ubiquityThreshold) + '-' + str(singleCopyThreshold) + '.tsv'
        fout = open(filename, 'w')
        fout.write('Genome ID\tLineage\tNumber of Genes\tUniformity\tDistribution\n')
        matrix = []
        rowLabels = []
        for genomeId in genomeIds:
            binSize = float(metadata[genomeId]['genome size']) / numBins

            binCounts = [0]*numBins
            pts = []
            for _, data in list(geneDistTable[genomeId].items()):
                for genePos in data:
                    binNum = int(genePos[1] / binSize)
                    binCounts[binNum] += 1
                    pts.append(genePos[1])
            matrix.append(binCounts)

            u = markerSet.uniformity(metadata[genomeId]['genome size'], pts)

            fout.write(genomeId + '\t' + '; '.join(metadata[genomeId]['taxonomy']) + '\t' + str(len(geneDistTable[genomeId])) + '\t%.3f' % u)
            for b in range(0, numBins):
                fout.write('\t' + str(binCounts[b]))
            fout.write('\n')

            rowLabels.append('%.2f' % u + ', ' + str(genomeId) + ' - ' + '; '.join(metadata[genomeId]['taxonomy'][0:5]))

        fout.close()

        # plot data
        heatmap = Heatmap()
        plotFilename = 'geneDistribution.' + taxonomyStr.replace(';','_') + '.' + str(ubiquityThreshold) + '-' + str(singleCopyThreshold) + '.png'
        heatmap.plot(plotFilename, matrix, rowLabels, 0.6)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate size of marker set.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-T', '--taxonomy', help='IMG taxonomy string indicating lineage of interest', default = 'universal')
    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-b', '--bins', help='Number of bins to divide genome into', type=int, default = 100)
    parser.add_argument('-n', '--num_genomes', help='Number of genomes in lineage to plot (-1 for all)', type=int, default = -1)

    args = parser.parse_args()

    plotMarkerSetDistribution = PlotMarkerSetDistribution()
    plotMarkerSetDistribution.run(args.taxonomy, args.ubiquity, args.single_copy, args.bins, args.num_genomes)
