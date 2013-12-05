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
Calculate marker set for varying lineages, number of genomes, and marker set parameters.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse
import random

from numpy import mean, std, arange

from lib.img import IMG
from lib.markerSet import MarkerGenes
from lib.plots.errorbar import ErrorBar
from lib.plots.boxplot import BoxPlot

class MarkerSetTest(object):
  def __init__(self):
    pass

  def run(self, taxonomyStr, ubiquityThreshold, singleCopyThreshold, replicates, minGenomes, maxGenomes, stepSize):
    img = IMG()
    markergenes = MarkerGenes()
    
    genomeIds = img.genomeIdsByTaxonomy(taxonomyStr, 'Final')
      
    print 'Lineage ' + taxonomyStr + ' contains ' + str(len(genomeIds)) + ' genomes.'
    if len(genomeIds) < minGenomes:
      sys.stderr.write('[Error] Insufficent number of genomes.\n')
      sys.exit()
      
    print ''
    print 'Ubiquity threshold: ' + str(ubiquityThreshold)
    print 'Single-copy threshold: ' + str(singleCopyThreshold)
      
    meanMarkerSetSize = []
    stdMarkerSetSize = []
    markerSetSizes = []
    if maxGenomes == -1:
      maxGenomes = len(genomeIds)
      
    if maxGenomes > len(genomeIds):
      maxGenomes = len(genomeIds)
      
    countTable = img.countTable(genomeIds) 
    countTable = img.filterTable(genomeIds, countTable)
   
    for numGenomes in xrange(minGenomes, maxGenomes, stepSize):
      markerSetSize = []
      for _ in xrange(0, replicates):
        genomeIdSubset = random.sample(genomeIds, numGenomes)
        
        markerGenes = markergenes.identify(genomeIdSubset, countTable, ubiquityThreshold*len(genomeIdSubset), singleCopyThreshold*len(genomeIdSubset))
        geneDistTable = img.geneDistTable(genomeIdSubset, markerGenes)
        colocatedGenes = img.colocatedGenes(geneDistTable)
        colocatedSets = img.colocatedSets(colocatedGenes, markerGenes)
        
        markerSetSize.append(len(colocatedSets))
        
      markerSetSizes.append(markerSetSize)
        
      m = mean(markerSetSize)
      meanMarkerSetSize.append(m)
      
      s = std(markerSetSize)
      stdMarkerSetSize.append(s)
        
      print ''
      print 'Genomes: ' + str(numGenomes) + ', Ubiquity > ' + str(int(ubiquityThreshold*len(genomeIdSubset))) + ', Single-copy > ' + str(int(singleCopyThreshold*len(genomeIdSubset)))
      print 'Mean: %.2f +/- %.2f' % (m, s)
      print 'Min: %d, Max: %d' %(min(markerSetSize), max(markerSetSize))
            
    # plot data
    errorBar = ErrorBar()
    plotFilename = './images/markerset.' + taxonomyStr.replace(';','_') + '.' + str(ubiquityThreshold) + '-' + str(singleCopyThreshold) +  '.errorbar.png'
    title = taxonomyStr.replace(';', '; ') + '\n' + 'Ubiquity = %.2f' % ubiquityThreshold + ', Single-copy = %.2f' % singleCopyThreshold
    errorBar.plot(plotFilename, arange(minGenomes, maxGenomes, stepSize), meanMarkerSetSize, stdMarkerSetSize, 'Number of Genomes', 'Marker Set Size', title)
    
    boxPlot = BoxPlot()
    plotFilename = './images/markerset.' + taxonomyStr.replace(';','_') + '.' + str(ubiquityThreshold) + '-' + str(singleCopyThreshold) +  '.boxplot.png'
    boxPlot.plot(plotFilename, markerSetSizes, arange(minGenomes, maxGenomes, stepSize), 'Number of Genomes', 'Marker Set Size', True, title)
    
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Plot change in marker set size for varying subsets of genomes.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-T', '--taxonomy', help='IMG taxonomy string indicating lineage of interest', default = 'universal')
  parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
  parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
  parser.add_argument('-r', '--replicates', help='Resampling replicates to perform for varying genomes sets', type=int, default=100)
  parser.add_argument('-i', '--min_genomes', help='Minimum genomes to consider when building marker set', type=int, default=10)
  parser.add_argument('-a', '--max_genomes', help='Maximum genomes to consider when building marker set, use -1 for all genomes', type=int, default=200)
  parser.add_argument('-x', '--step_size', help='Step size for increasing number of genomes under consideration', type=int, default=5)

  args = parser.parse_args()
  
  markerSetTest = MarkerSetTest()
  markerSetTest.run(args.taxonomy, args.ubiquity, args.single_copy, args.replicates, args.min_genomes, args.max_genomes, args.step_size)

