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
Plot scaffold length vs. percentage of marker genes.
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
from collections import defaultdict

from checkm.taxonParser import TaxonParser
from checkm.util.img import IMG
from checkm.plot.AbstractPlot import AbstractPlot

class PlotOptions:
    def __init__(self):
        self.font_size = 8
        self.dpi = 300
        self.width = 6
        self.height = 6

class ScatterPlot(AbstractPlot):
    def __init__(self):
        plotOptions = PlotOptions()
        AbstractPlot.__init__(self, plotOptions)
    
    def plot(self, x, y): 
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axes = self.fig.add_subplot(111)

        self.plotOnAxes(axes, x, y)
        
        #self.fig.tight_layout()
        self.draw()
        
    def plotOnAxes(self, axes, x, y):
        axes.scatter(x, y, s=10, lw=0.5)    
        axes.set_xlabel('Scaffold length (kbp)')
        axes.set_ylabel('% marker genes')
                
        # Change sequence lengths from bp to kbp
        axes.set_xscale('log')
        xticks = axes.get_xticks()
        kbpLabels = []
        for seqLen in xticks:
            label = '%.1f' % (float(seqLen)/1000)
            label = label.replace('.0', '') # remove trailing zero
            kbpLabels.append(label)
        axes.set_xticklabels(kbpLabels)
            
        # Prettify plot     
        for a in axes.yaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
                
        for a in axes.xaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
            
        for line in axes.yaxis.get_ticklines(): 
            line.set_color(self.axesColour)
                
        for line in axes.xaxis.get_ticklines(): 
            line.set_color(self.axesColour)
            
        for loc, spine in axes.spines.iteritems():
            if loc in ['right','top']:
                spine.set_color('none') 
            else:
                spine.set_color(self.axesColour)


class PlotScaffoldLenVsMarkers(object):
    def __init__(self):
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        

    def run(self):
        # get all draft genomes consisting of a user-specific minimum number of scaffolds
        print ''
        metadata = self.img.genomeMetadata()
        print '  Total genomes: %d' % len(metadata)
        
        arGenome = set()
        for genomeId in metadata:
            if metadata[genomeId]['taxonomy'][0] == 'Archaea':
                arGenome.add(genomeId)
                
        draftGenomeIds = arGenome - self.img.filterGenomeIds(arGenome, metadata, 'status', 'Finished')
        print '  Number of draft genomes: %d' % len(draftGenomeIds)
        
        minScaffolds = 20
        genomeIdsToTest = set()
        for genomeId in draftGenomeIds:
            if metadata[genomeId]['scaffold count'] >= minScaffolds:
                genomeIdsToTest.add(genomeId)
        print '  Number of draft genomes with >= %d scaffolds: %d' % (minScaffolds, len(genomeIdsToTest))

        print ''
        print '  Calculating genome information for calculating marker sets:'
        genomeFamilyScaffolds = self.img.precomputeGenomeFamilyScaffolds(genomeIdsToTest)
        
        print '  Calculating genome sequence lengths.'
        genomeSeqLens = self.img.precomputeGenomeSeqLens(genomeIdsToTest)
        
        print '  Determining domain-specific marker sets.'
        taxonParser = TaxonParser()
        taxonMarkerSets = taxonParser.readMarkerSets()
        bacMarkers = taxonMarkerSets['domain']['Bacteria'].getMarkerGenes()
        arMarkers = taxonMarkerSets['domain']['Archaea'].getMarkerGenes()
        print '    There are %d bacterial markers and %d archaeal markers.' % (len(bacMarkers), len(arMarkers))
        
        print '  Determining percentage of markers on each scaffold.'
        totalMarkers = 0
        totalSequenceLen = 0
        markersOnShortScaffolds = 0
        totalShortScaffoldLen = 0
        
        scaffoldLen = {}
        percentageMarkers = defaultdict(float)
        for genomeId, markerIds in genomeFamilyScaffolds.iteritems():
            domain = metadata[genomeId]['taxonomy'][0]
            markerGenes = bacMarkers if domain == 'Bacteria' else arMarkers
            for markerId in markerGenes:
                if markerId.startswith('PF'):
                    markerId = markerId.replace('PF', 'pfam')
                    markerId = markerId[0:markerId.rfind('.')]
                if markerId in markerIds:
                    for scaffoldId in markerIds[markerId]:
                        scaffoldLen[scaffoldId] = genomeSeqLens[genomeId][scaffoldId]
                        percentageMarkers[scaffoldId] += 1.0/len(markerGenes)
                        
                        totalMarkers += 1
                        totalSequenceLen += genomeSeqLens[genomeId][scaffoldId]
                        
                        if genomeSeqLens[genomeId][scaffoldId] < 10000:
                            markersOnShortScaffolds += 1
                            totalShortScaffoldLen += genomeSeqLens[genomeId][scaffoldId]
       
        print 'Markers on short scaffolds: %d over %d Mbp (%f markers per base)' % (markersOnShortScaffolds, totalShortScaffoldLen, float(markersOnShortScaffolds)/totalShortScaffoldLen)
        print 'Total markers on scaffolds: %d over %d Mbp (%f markers per base)' % (totalMarkers, totalSequenceLen, float(totalMarkers)/totalSequenceLen)
                        
        print '  Create plot.'
        plotLens = []
        plotPerMarkers = []
        for scaffoldId in percentageMarkers:
            plotLens.append(scaffoldLen[scaffoldId])
            plotPerMarkers.append(percentageMarkers[scaffoldId]/scaffoldLen[scaffoldId] * 1e6)
            
        scatterPlot = ScatterPlot()
        scatterPlot.plot(plotLens, plotPerMarkers)     
        scatterPlot.savePlot('./experiments/plotScaffoldLenVsMarkers.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    plotScaffoldLenVsMarkers = PlotScaffoldLenVsMarkers()
    plotScaffoldLenVsMarkers.run()
