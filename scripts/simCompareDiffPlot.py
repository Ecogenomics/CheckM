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
Plot showing the percentage of times the lineage-specific marker set
  outperforms the domain-specific marker set when they deviate from 
  each other by different amounts.
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

from checkm.util.img import IMG
from checkm.plot.AbstractPlot import AbstractPlot

import numpy as np

class PlotOptions:
    def __init__(self):
        self.font_size = 8
        self.dpi = 1200
        self.width = 6
        self.height = 4

class StackedBarPlot(AbstractPlot):
    def __init__(self):
        plotOptions = PlotOptions()
        AbstractPlot.__init__(self, plotOptions)
    
    def plot(self, lineageCountsComp, domainCountsComp, lineageCountsCont, domainCountsCont): 
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesComp = self.fig.add_subplot(211)
        self.plotOnAxes(axesComp, lineageCountsComp, domainCountsComp, False, True)
        
        axesCont = self.fig.add_subplot(212)
        self.plotOnAxes(axesCont, lineageCountsCont, domainCountsCont, True, False)
        
        self.fig.tight_layout(pad=5, h_pad=15)
        self.draw()
        
    def plotOnAxes(self, axes, lineageCounts, domainCounts, bXLabel, bLegend):  
        width = 0.80
        index = np.arange(len(lineageCounts))
        p1 = axes.bar(index, lineageCounts, width, color='0.4', lw=0)  
        p2 = axes.bar(index, domainCounts, width, color='0.7', bottom=lineageCounts, lw=0)    
        
        if bXLabel:
            axes.set_xlabel('% differences between domain and lineage estimates')
        axes.set_ylabel('% test genomes')
                
        labels = []
        for i in xrange(0, len(lineageCounts), 2):
            labels.append('%d%%' % (i+1))
            labels.append('')
                  
        axes.set_xticks(index+width/2)
        axes.set_xticklabels(labels)
        axes.set_xlim(0, len(lineageCounts))
        
        if bLegend:
            legend = axes.legend( (p2[0], p1[0]), ('Domain', 'Lineage') )
            legend.draw_frame(False)
          
        # Prettify plot     
        for a in axes.yaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
                
        for a in axes.xaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
            
        for line in axes.yaxis.get_ticklines(): 
            line.set_color('black')
                
        for line in axes.xaxis.get_ticklines(): 
            line.set_color('black')
            
        for loc, spine in axes.spines.iteritems():
            if loc in ['right','top']:
                spine.set_color('none') 
            else:
                spine.set_color('black')

class SimCompareDiffPlot(object):
    def __init__(self):
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        
    def run(self):
        # count number of times the lineage-specific marker set results outperform
        # the domain-specific marker set for varying differences between the two sets
        numBars = 15
        
        lineageCountsComp = [0]*numBars
        domainCountsComp = [0]*numBars
        
        lineageCountsCont = [0]*numBars
        domainCountsCont = [0]*numBars
        
        totalCountsComp = 0
        totalCountsCont = 0
        
        domCompBest = 0
        lineageCompBest = 0
        domContBest = 0
        lineageContBest = 0
        
        metadata = self.img.genomeMetadata()
        domCompTaxon = defaultdict(int)
        lineageCompTaxon = defaultdict(int)
        
        for line in open('./simulations/briefSummaryOut.tsv'):
            lineSplit = line.split('\t')
            genomeId = lineSplit[0]
            taxonomy = metadata[genomeId]['taxonomy']
            phylum = taxonomy[1]
            domCompMS, lineageCompMS, lineageCompRMS, domContMS, lineageContMS, lineageContRMS = [float(x) for x in lineSplit[1:]]
            
            diff = abs(abs(lineageCompMS) - abs(domCompMS))
            if diff > 5:
                intDiff = int(diff)
                if intDiff >= numBars:
                    intDiff = (numBars-1)
                    
                if abs(domCompMS) < abs(lineageCompMS):
                    domainCountsComp[intDiff] += 1
                    domCompBest += 1
                    domCompTaxon[phylum] += 1
                else:
                    lineageCountsComp[intDiff] += 1
                    lineageCompBest += 1
                    lineageCompTaxon[phylum] += 1
                    
                totalCountsComp += 1
                
            diff = abs(abs(lineageContMS) - abs(domContMS))
            if diff > 5:
                intDiff = int(diff)
                if intDiff >= numBars:
                    intDiff = (numBars-1)
                    
                if abs(domContMS) < abs(lineageContMS):
                    domainCountsCont[intDiff] += 1
                    domContBest += 1
                else:
                    lineageCountsCont[intDiff] += 1
                    lineageContBest += 1
                    
                totalCountsCont += 1
                
        print '%% times lineage comp better than domain: %.2f' % (float(lineageCompBest)*100/(domCompBest + lineageCompBest))
        print '%% times lineage cont better than domain: %.2f' % (float(lineageContBest)*100/(domContBest + lineageContBest))
        
        print ''
        print 'Taxonomy breakdown (dom best, lineage best):'
        taxa = set(domCompTaxon.keys()).union(lineageCompTaxon.keys())
        for t in taxa:
            print '%s\t%.2f\t%.2f' % (t, domCompTaxon[t]*100.0/domCompBest, lineageCompTaxon[t]*100.0/lineageCompBest)
                
        # normalize counts
        for i in xrange(0, numBars):
            lineageCountsComp[i] = float(lineageCountsComp[i])*100 / totalCountsComp
            domainCountsComp[i] = float(domainCountsComp[i])*100 / totalCountsComp
            
            if domainCountsComp[i] > lineageCountsComp[i]:
                print 'Domain bets lineage (comp): %d%% (%f, %f)' % (i+1, domainCountsComp[i], lineageCountsComp[i])
            
            lineageCountsCont[i] = float(lineageCountsCont[i])*100 / totalCountsCont
            domainCountsCont[i] = float(domainCountsCont[i])*100 / totalCountsCont
            
            if domainCountsCont[i] > lineageCountsCont[i]:
                print 'Domain bets lineage (cont): %d%% (%f, %f)' % (i+1, domainCountsCont[i], lineageCountsCont[i])
         
        stackedBarPlot = StackedBarPlot()
        stackedBarPlot.plot(lineageCountsComp, domainCountsComp, lineageCountsCont, domainCountsCont)     
        stackedBarPlot.savePlot('./experiments/simCompareDiffPlot.svg')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simCompareDiffPlot = SimCompareDiffPlot()
    simCompareDiffPlot.run()
