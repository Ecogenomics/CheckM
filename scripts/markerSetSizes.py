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
Calculate table of marker set for all lineages.
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

from numpy import arange

from lib.img import IMG

class MarkerSetSizes(object):
    def __init__(self):
        pass

    def run(self, minThreshold, maxThreshold, stepSize, minGenomes, mostSpecificRanks):
        img = IMG()

        trustedGenomeIds = img.trustedGenomes()

        fout = open('./data/markerSetSize.tsv', 'w')
        fout.write('Lineage\t# genomes')
        for threshold in arange(maxThreshold, minThreshold, -stepSize):
            fout.write('\t' + str(threshold))
        fout.write('\n')

        lineages = img.lineagesSorted(mostSpecificRanks)
        for lineage in lineages:
            genomeIds = img.genomeIdsByTaxonomy(lineage)
            genomeIds = list(genomeIds.intersection(trustedGenomeIds))

            if len(genomeIds) < minGenomes:
                continue

            print(('\nLineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes.'))
            fout.write(lineage + '\t' + str(len(genomeIds)))

            pfamTable = img.pfamTable(genomeIds)
            for threshold in arange(maxThreshold, minThreshold, -stepSize):
                markerSet = img.markerGenes(genomeIds, pfamTable, threshold*len(genomeIds), threshold*len(genomeIds))
                fout.write('\t' + str(len(markerSet)))
                print(('  Threshold = %.2f, marker set size = %d' % (threshold, len(markerSet))))
            fout.write('\n')

        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate table of marker set for all lineage.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a', '--min_threshold', help='Minimum threshold for defining marker set', type=float, default = 0.90)
    parser.add_argument('-b', '--max_threshold', help='Maximum threshold for defining marker set', type=float, default = 1.0)
    parser.add_argument('-x', '--step_size', help='Step size for increasing number of genomes under consideration', type=int, default=0.01)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 10)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 5)

    args = parser.parse_args()

    markerSetSizes = MarkerSetSizes()
    markerSetSizes.run(args.min_threshold, args.max_threshold, args.step_size, args.min_genomes, args.most_specific_rank)
