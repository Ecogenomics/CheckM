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
Calculate size of marker set for specific lineages.
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
from lib.plots.lineplot import LinePlot

class MarkerSetTest(object):
  def __init__(self):
    pass

  def run(self, taxonomyStr, minThreshold, maxThreshold, stepSize):
    img = IMG()
    
    genomeIds = img.genomeIdsByTaxonomy(taxonomyStr, 'Final')
      
    print 'Lineage ' + taxonomyStr + ' contains ' + str(len(genomeIds)) + ' genomes.'
      
    markerSetSizes = []
   
    countTable = img.countTable(genomeIds)
    for threshold in arange(maxThreshold, minThreshold, -stepSize):
      markerGenes = img.markerGenes(genomeIds, countTable, threshold*len(genomeIds), threshold*len(genomeIds))
      
      geneDistTable = img.geneDistTable(genomeIds, markerGenes)
      colocatedGenes = img.colocatedGenes(geneDistTable)
      colocatedSets = img.colocatedSets(colocatedGenes, markerGenes)
      
      markerSetSizes.append(len(colocatedSets))
      
      print '  Threshold = %.2f, marker set size = %d' % (threshold, len(markerGenes))

    # plot data
    plot = LinePlot()
    plotFilename = './images/markerSetSize.' + taxonomyStr.replace(';','_') + '.png'
    title = taxonomyStr.replace(';', '; ')
    plot.plot(plotFilename, arange(maxThreshold, minThreshold, -stepSize), markerSetSizes, 'Threshold', 'Marker Set Size', title)
    
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Calculate size of marker set for specific lineage.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-T', '--taxonomy', help='IMG taxonomy string indicating lineage of interest', default = 'prokaryotes')
  parser.add_argument('-a', '--min_threshold', help='Minimum threshold for defining marker set', type=float, default = 0.85)
  parser.add_argument('-b', '--max_threshold', help='Maximum threshold for defining marker set', type=float, default = 1.0)
  parser.add_argument('-x', '--step_size', help='Step size for increasing number of genomes under consideration', type=int, default=0.01)

  args = parser.parse_args()
  
  markerSetTest = MarkerSetTest()
  markerSetTest.run(args.taxonomy, args.min_threshold, args.max_threshold, args.step_size)

