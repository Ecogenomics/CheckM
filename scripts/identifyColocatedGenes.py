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
Identify genes that are consistently co-located within genomes from a specific lineage.
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

from lib.img import IMG
from lib.markerSet import MarkerSet

class IdentifyColocatedGenes(object):
    def __init__(self):
        pass

    def run(self, ubiquityThreshold, singleCopyThreshold, minGenomes, minMarkers, mostSpecificRank, distThreshold, genomeThreshold):
        img = IMG()
        markerset = MarkerSet()

        lineages = img.lineagesSorted(mostSpecificRank)

        fout = open('./data/colocated.tsv', 'w', 1)
        fout.write('Lineage\t# genomes\t# markers\t# co-located sets\tCo-located markers\n')

        lineageCount = 0
        for lineage in lineages:
            lineageCount += 1

            genomeIds = img.genomeIdsByTaxonomy(lineage, 'Final')
            if len(genomeIds) < minGenomes:
                continue

            countTable = img.countTable(genomeIds)
            markerGenes = markerset.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))

            geneDistTable = img.geneDistTable(genomeIds, markerGenes)
            colocatedGenes = markerset.colocatedGenes(geneDistTable, distThreshold, genomeThreshold)
            colocatedSets = markerset.colocatedSets(colocatedGenes, markerGenes)
            if len(colocatedSets) < minMarkers:
                continue

            print '\nLineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes (' + str(lineageCount) + ' of ' + str(len(lineages)) + ').'
            print '  Marker genes: ' + str(len(markerGenes))
            print '  Co-located gene sets: ' + str(len(colocatedSets))

            fout.write(lineage + '\t' + str(len(genomeIds)) + '\t' + str(len(markerGenes)) + '\t' + str(len(colocatedSets)))
            for cs in colocatedSets:
                fout.write('\t' + ', '.join(cs))
            fout.write('\n')

        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify co-located genes within genomes for a specific lineage.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 40)
    parser.add_argument('-x', '--min_markers', help='Minimum markers required to include in analysis', type=int, default = 20)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 5)
    parser.add_argument('-d', '--distance_threshold', help='Distance, as percentage of genome length, to be considered co-located', type=float, default=0.001)
    parser.add_argument('-g', '--genome_threshold', help='Percentage of genomes required to be considered co-located', type=float, default=0.95)

    args = parser.parse_args()

    identifyColocatedGenes = IdentifyColocatedGenes()
    identifyColocatedGenes.run(args.ubiquity, args.single_copy, args.min_genomes, args.min_markers, args.most_specific_rank, args.distance_threshold, args.genome_threshold)
