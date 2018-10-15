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
Perform simulation to show that marker sets give better completion estimations
  compared to marker genes.
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
from lib.taxonomyUtils import ranksByLabel
from lib.plots.boxplot import BoxPlot

class SimMarkerGenesVsMarkerSets(object):
    def __init__(self):
        pass

    def run(self, taxonomyStr, ubiquityThreshold, singleCopyThreshold, percentCompletion, numReplicates, numGenomes, contigLen):
        img = IMG()

        genomeIds = img.genomeIdsByTaxonomy(taxonomyStr, 'Final')
        print(('\nLineage ' + taxonomyStr + ' contains ' + str(len(genomeIds)) + ' genomes.'))

        # build marker genes and colocated marker sets
        countTable = img.countTable(genomeIds)
        markerGenes = img.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
        print(('  Marker genes: ' + str(len(markerGenes))))

        geneDistTable = img.geneDistTable(genomeIds, markerGenes, spacingBetweenContigs=1e6)
        colocatedGenes = img.colocatedGenes(geneDistTable)
        colocatedSets = img.colocatedSets(colocatedGenes, markerGenes)
        print(('  Co-located gene sets: ' + str(len(colocatedSets))))


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
            mgCompletion = []
            msCompletion = []
            for _ in range(0, numReplicates):
                startPartialGenomeContigs = img.sampleGenome(metadata[genomeId]['genome size'], percentCompletion, contigLen)

                # calculate completion with marker genes
                containedMarkerGenes = img.containedMarkerGenes(markerGenes, geneDistTable[genomeId], startPartialGenomeContigs, contigLen)
                mgCompletion.append(float(len(containedMarkerGenes))/len(markerGenes) - percentCompletion)

                # calculate completion with marker set
                comp = 0.0
                for cs in colocatedSets:
                    present = 0
                    for contigId in cs:
                        if contigId in containedMarkerGenes:
                            present += 1

                    comp += float(present) / len(cs)
                msCompletion.append(comp / len(colocatedSets) - percentCompletion)

            plotData.append(mgCompletion)
            plotData.append(msCompletion)

            species = ' '.join(metadata[genomeId]['taxonomy'][ranksByLabel['Genus']:])

            plotLabels.append(species + ' (' + genomeId + ')')
            plotLabels.append('')

        # plot data
        boxPlot = BoxPlot()
        plotFilename = './images/sim.MGvsMS.' + taxonomyStr.replace(';','_') + '.' + str(percentCompletion) + '.errorbar.png'
        title = taxonomyStr.replace(';', '; ') + '\n' + 'Percent completion = %.2f' % percentCompletion
        boxPlot.plot(plotFilename, plotData, plotLabels, r'$\Delta$' + ' Percent Completion', '', False, title)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-T', '--taxonomy', help='IMG taxonomy string indicating lineage of interest', default = 'prokaryotes')
    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-p', '--percent_complete', help='Percent completion to simulate', type=float, default = 0.75)
    parser.add_argument('-r', '--replicates', help='Replicates per genome.', type=int, default = 100)
    parser.add_argument('-g', '--num_genomes', help='Number of random genomes to consider (-1 for all)', type=int, default = 20)
    parser.add_argument('-c', '--contig_len', help='Length of contigs to simulate', type=int, default = 5000)

    args = parser.parse_args()

    simMarkerGenesVsMarkerSets = SimMarkerGenesVsMarkerSets()
    simMarkerGenesVsMarkerSets.run(args.taxonomy, args.ubiquity, args.single_copy, args.percent_complete, args.replicates, args.num_genomes, args.contig_len)
