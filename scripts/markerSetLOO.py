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
Calculate marker set for varying lineages under a leave-one-out scheme.
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

from numpy import mean, std

from lib.img import IMG
from lib.plots.boxplot import BoxPlot

class MarkerSetLOO(object):
    def __init__(self):
        pass

    def run(self, ubiquityThreshold, singleCopyThreshold, minGenomes, mostSpecificRank, minMarkers):
        print(('Ubiquity threshold: ' + str(ubiquityThreshold)))
        print(('Single-copy threshold: ' + str(singleCopyThreshold)))
        print(('Min. genomes: ' + str(minGenomes)))
        print(('Most specific taxonomic rank: ' + str(mostSpecificRank)))

        img = IMG()

        deltaMarkerSetSizes = []

        lineages = img.lineagesByCriteria(minGenomes, mostSpecificRank)
        lineages = ['prokaryotes'] + lineages

        boxPlotLabels = []
        for lineage in lineages:
            genomeIds = img.genomeIdsByTaxonomy(lineage)
            trusted = img.trustedGenomes()
            genomeIds = list(genomeIds.intersection(trusted))

            print('')
            print(('Lineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes.'))

            # get table of PFAMs and do some initial filtering to remove PFAMs that are
            # clearly not going to pass the ubiquity and single-copy thresholds
            pfamTable = img.pfamTable(genomeIds)
            pfamTable = img.filterPfamTable(genomeIds, pfamTable, ubiquityThreshold*0.9, singleCopyThreshold*0.9)

            markerSet = img.markerGenes(genomeIds, pfamTable, ubiquityThreshold*(len(genomeIds)-1), singleCopyThreshold*(len(genomeIds)-1))
            fullMarkerSetSize = len(markerSet)

            if fullMarkerSetSize < minMarkers:
                continue

            boxPlotLabels.append(lineage.split(';')[-1].strip() + ' (' + str(len(genomeIds)) + ', ' + str(fullMarkerSetSize) + ')')

            deltaMarkerSetSize = []
            numGenomes = len(genomeIds)-1

            for loo in range(0, len(genomeIds)):
                if loo != len(genomeIds) - 1:
                    genomeIdSubset = genomeIds[0:loo] + genomeIds[loo+1:]
                else:
                    genomeIdSubset = genomeIds[0:loo]

                markerSet = img.markerGenes(genomeIdSubset, pfamTable, ubiquityThreshold*len(genomeIdSubset), singleCopyThreshold*len(genomeIdSubset))
                deltaMarkerSetSize.append(fullMarkerSetSize - len(markerSet))

                if fullMarkerSetSize < len(markerSet):
                    print('[Warning] Unexpected!')

            deltaMarkerSetSizes.append(deltaMarkerSetSize)

            m = mean(deltaMarkerSetSize)
            s = std(deltaMarkerSetSize)

            print(('  LOO Ubiquity >= ' + str(int(ubiquityThreshold*numGenomes)) + ', LOO Single-copy >= ' + str(int(singleCopyThreshold*numGenomes))))
            print(('  Delta Mean: %.2f +/- %.2f' % (m, s)))
            print(('  Delta Min: %d, Delta Max: %d' % (min(deltaMarkerSetSize), max(deltaMarkerSetSize))))

        # plot data
        boxPlot = BoxPlot()
        plotFilename = './images/LOO.' + str(ubiquityThreshold) + '-' + str(singleCopyThreshold) + '.boxplot.png'
        title = 'Ubiquity = %.2f' % ubiquityThreshold + ', Single-copy = %.2f' % singleCopyThreshold
        boxPlot.plot(plotFilename, deltaMarkerSetSizes, boxPlotLabels, r'$\Delta$' + ' Marker Set Size', '', False, title)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate marker set for varying lineages under a leave-one-out scheme.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 20)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 5)
    parser.add_argument('-x', '--min_markers', help='Minimum markers required to include in analysis', type=int, default = 20)

    args = parser.parse_args()

    markerSetLOO = MarkerSetLOO()
    markerSetLOO.run(args.ubiquity, args.single_copy, args.min_genomes, args.most_specific_rank, args.min_markers)
