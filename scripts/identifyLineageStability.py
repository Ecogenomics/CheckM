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
Identify lineages with sufficiently stable marker sets to use for completion and
  contamination estimation.
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

from numpy import mean, std

class IdentifySufficientLineages(object):
    def __init__(self):
        pass

    def run(self, ubiquityThreshold, singleCopyThreshold, minGenomes, minMarkers, mostSpecificRank, percentGenomes, numReplicates):
        img = IMG()

        lineages = img.lineagesByCriteria(minGenomes, mostSpecificRank)

        fout = open('./data/lineage_evaluation.tsv', 'w')
        fout.write('Lineage\t# genomes\t# markers\tpercentage\tnum replicates\tmean\tstd\tmean %\tmean + std%\tmean + 2*std %\n')

        for lineage in lineages:
            genomeIds = img.genomeIdsByTaxonomy(lineage, 'Final')
            if len(genomeIds) < minGenomes:
                continue

            countTable = img.countTable(genomeIds)
            countTable = img.filterTable(genomeIds, countTable, ubiquityThreshold*0.9, singleCopyThreshold*0.9)

            # calculate marker set for all genomes
            markerGenes = img.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
            if len(markerGenes) < minMarkers:
                continue

            print(('\nLineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes.'))
            print(('  Marker genes: ' + str(len(markerGenes))))

            fout.write(lineage + '\t' + str(len(genomeIds)) + '\t' + str(len(markerGenes)) + '\t%.2f' % percentGenomes + '\t' + str(numReplicates))

            # withhold select percentage of genomes and calculate new marker set
            changeMarkerSetSize = []
            for _ in range(0, numReplicates):
                subsetGenomeIds = random.sample(genomeIds, int((1.0-percentGenomes)*len(genomeIds) + 0.5))

                newMarkerGenes = img.markerGenes(subsetGenomeIds, countTable, ubiquityThreshold*len(subsetGenomeIds), singleCopyThreshold*len(subsetGenomeIds))

                changeMarkerSetSize.append(len(newMarkerGenes.symmetric_difference(markerGenes)))

            m = mean(changeMarkerSetSize)
            s = std(changeMarkerSetSize)

            print(('  Mean: %.2f, Std: %.2f, Per: %.2f' % (m, s, (m+ 2*s) * 100 / len(markerGenes))))
            fout.write('\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (m, s, m * 100 / len(markerGenes), (m + s) * 100 / len(markerGenes), (m + 2*s) * 100 / len(markerGenes)) + '\n')

        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify co-located genes within genomes for a specific lineage.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 20)
    parser.add_argument('-x', '--min_markers', help='Minimum markers required to include in analysis', type=int, default = 20)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 5)

    parser.add_argument('-p', '--percent_genomes', help='Percentage of genomes to withhold', type=float, default = 0.1)
    parser.add_argument('-y', '--replicates', help='Replicates to perform', type=int, default = 100)

    args = parser.parse_args()

    identifySufficientLineages = IdentifySufficientLineages()
    identifySufficientLineages.run(args.ubiquity, args.single_copy, args.min_genomes, args.min_markers, args.most_specific_rank, args.percent_genomes, args.replicates)
