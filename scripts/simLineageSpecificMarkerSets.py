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
Perform simulation to show that lineage specific marker sets give better completion
  estimations compared to more general marker sets.
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

from lib.img import IMG
from lib.plots.boxplot import BoxPlot

class SimLineageSpecificMarkerSets(object):
    def __init__(self):
        pass

    def run(self, taxonomyStr, mostSpecificRank, minGenomes, ubiquityThreshold, singleCopyThreshold, percentCompletion, numReplicates, numGenomes, contigLen):
        img = IMG()

        lineages = []
        taxon = taxonomyStr.split(';')
        for r in range(0, len(taxon)):
            lineages.append(';'.join(taxon[0:r+1]))

        # get all marker sets
        markerGenes = []
        geneDistTable = []
        colocatedSets = []
        for lineage in lineages:
            genomeIds = img.genomeIdsByTaxonomy(lineage, 'Final')
            print('\nLineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes.')

            # build marker genes and colocated marker sets
            countTable = img.countTable(genomeIds)
            mg = img.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
            print('  Marker genes: ' + str(len(mg)))

            mdt = img.geneDistTable(genomeIds, mg, spacingBetweenContigs=1e6)
            colocatedGenes = img.colocatedGenes(mdt)
            cs = img.colocatedSets(colocatedGenes, mg)
            print('  Co-located gene sets: ' + str(len(cs)))

            markerGenes.append(mg)
            geneDistTable.append(mdt)
            colocatedSets.append(cs)

        # random sample genomes
        if numGenomes == -1:
            rndGenomeIds = genomeIds
        else:
            rndGenomeIds = random.sample(genomeIds, numGenomes)

        # estimate completion for each genome using both the marker genes and marker sets
        metadata = img.genomeMetadata('Final')
        plotLabels = []
        plotData = []
        for genomeId in rndGenomeIds:
            completion = [[] for _ in range(len(lineages))]
            for _ in range(0, numReplicates):
                startPartialGenomeContigs = img.sampleGenome(metadata[genomeId]['genome size'], percentCompletion, contigLen)

                # calculate completion with marker set
                for i in range(len(lineages)):
                    containedMarkerGenes = img.containedMarkerGenes(markerGenes[i], geneDistTable[i][genomeId], startPartialGenomeContigs, contigLen)

                    comp = 0.0
                    for cs in colocatedSets[i]:
                        present = 0
                        for contigId in cs:
                            if contigId in containedMarkerGenes:
                                present += 1

                        comp += float(present) / len(cs)

                    completion[i].append(comp / len(colocatedSets[i]) - percentCompletion)

                    plotLabels.append(genomeId + '  - ' + lineages[i])

            for d in completion:
                plotData.append(d)

        # plot data
        boxPlot = BoxPlot()
        plotFilename = './images/sim.lineages.' + taxonomyStr.replace(';','_') + '.' + str(percentCompletion) + '.errorbar.png'
        title = taxonomyStr.replace(';', '; ') + '\n' + 'Percent completion = %.2f' % percentCompletion
        boxPlot.plot(plotFilename, plotData, plotLabels, r'$\Delta$' + ' Percent Completion', '', False, title)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-T', '--taxonomy', help='IMG taxonomy string indicating lineage of interest', default = 'Archaea; Euryarchaeota')
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 1)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 20)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-p', '--percent_complete', help='Percent completion to simulate', type=float, default = 0.75)
    parser.add_argument('-x', '--replicates', help='Replicates per genome.', type=int, default = 100)
    parser.add_argument('-g', '--num_genomes', help='Number of random genomes to consider (-1 for all)', type=int, default = 20)
    parser.add_argument('-c', '--contig_len', help='Length of contigs to simulate', type=int, default = 5000)

    args = parser.parse_args()

    simLineageSpecificMarkerSets = SimLineageSpecificMarkerSets()
    simLineageSpecificMarkerSets.run(args.taxonomy, args.most_specific_rank, args.min_genomes, args.ubiquity, args.single_copy, args.percent_complete, args.replicates, args.num_genomes, args.contig_len)
