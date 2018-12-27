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
Calculate conservative list of marker genes for each lineage.
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

class CalculateConservativeMarkerGenes(object):
    def __init__(self):
        pass

    def run(self, ubiquityThreshold, singleCopyThreshold, rank):
        img = IMG()
        markerset = MarkerSet()

        print('Reading metadata.')
        metadata = img.genomeMetadata()
        print('  Genomes with metadata: ' + str(len(metadata)))

        # calculate marker set for each lineage at the specified rank
        sortedLineages = img.lineagesSorted(metadata, rank)
        markerGeneLists = {}
        for lineage in sortedLineages:
            taxonomy = lineage.split(';')
            if len(taxonomy) != rank+1:
                continue

        genomeIds = img.genomeIdsByTaxonomy(lineage, metadata, 'Final')
        countTable = img.countTable(genomeIds)

        if len(genomeIds) < 3:
            continue

        print('Lineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes.')

        markerGenes = markerset.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))

        print('  Marker genes: ' + str(len(markerGenes)))
        print('')

        markerGeneLists[lineage] = markerGenes

        # calculate union of marker gene list for higher taxonomic groups
        for r in range(rank-1, -1, -1):
            print('Processing rank ' + str(r))
            rankMarkerGeneLists = {}
            for lineage, markerGenes in markerGeneLists.items():
                taxonomy = lineage.split(';')
                if len(taxonomy) != r+2:
                    continue

                curLineage = '; '.join(taxonomy[0:r+1])
                if curLineage not in rankMarkerGeneLists:
                    rankMarkerGeneLists[curLineage] = markerGenes
                else:
                    curMarkerGenes = rankMarkerGeneLists[curLineage]
                    curMarkerGenes = curMarkerGenes.intersection(markerGenes)
                    rankMarkerGeneLists[curLineage] = curMarkerGenes

            # combine marker gene list dictionaries
            markerGeneLists.update(rankMarkerGeneLists)

    print('Archaeal markers: ' + str(len(markerGeneLists['Archaea'])))
    print('Bacterial markers: ' + str(len(markerGeneLists['Bacteria'])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate conservative list of marker genes for each lineage.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-r', '--rank', help='Rank at which to build initial marker genes', type=int, default = 3)

    args = parser.parse_args()

    calculateConservativeMarkerGenes = CalculateConservativeMarkerGenes()
    calculateConservativeMarkerGenes.run(args.ubiquity, args.single_copy, args.rank)
