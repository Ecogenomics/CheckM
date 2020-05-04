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
Histograms showing bias in estimations.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

from checkm.plot.AbstractPlot import AbstractPlot

import numpy as np

import matplotlib.mlab as mlab

class PlotOptions:
    def __init__(self):
        self.font_size = 8
        self.dpi = 300
        self.width = 8
        self.height = 8

class Histogram(AbstractPlot):
    def __init__(self):
        plotOptions = PlotOptions()
        AbstractPlot.__init__(self, plotOptions)
    
    def plot(self): 
        # read data
        fin = open('./simulations/correct_comp_cont_0.70-0.20.cont_per_by_g2.comp_by_im.tsv')
        data = fin.readlines()
        fin.close()
        
        
        compDom = map(float, data[0].split('\t'))
        contDom =map(float, data[1].split('\t'))
        compLineage = map(float, data[2].split('\t'))
        contLineage =map(float, data[3].split('\t'))
        
        correctCompDom = map(float, data[4].split('\t'))
        correctContDom =map(float, data[5].split('\t'))
        correctCompLineage = map(float, data[6].split('\t'))
        correctContLineage =map(float, data[7].split('\t'))
        
        compMin = min(compDom + compLineage + correctCompDom + correctCompLineage)
        compMax = max(compDom + compLineage + correctCompDom + correctCompLineage)
        
        contMin = min(contDom + contLineage + correctContDom + correctContLineage)
        contMax = max(contDom + contLineage + correctContDom + correctContLineage)
                    
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesComp = self.fig.add_subplot(221)
        self.plotOnAxes(axesComp, compDom, correctCompDom, compMin, compMax, 0, 0.1, '% error (comp domain)')
        
        axesCont = self.fig.add_subplot(222)
        self.plotOnAxes(axesCont, contDom, correctContDom, contMin, contMax, 0, 0.1, '% error (cont domain)')
        
        axesComp = self.fig.add_subplot(223)
        self.plotOnAxes(axesComp, compLineage, correctCompLineage, compMin, compMax, 0, 0.1, '% error (comp lineage)')
        
        axesCont = self.fig.add_subplot(224)
        self.plotOnAxes(axesCont, contLineage, correctContLineage, contMin, contMax, 0, 0.1, '% error (cont lineage)')
        
        self.fig.tight_layout(pad=2, h_pad=5)
        self.draw()
        
    def plotOnAxes(self, axes, beforeCorrection, afterCorrection, xMin, xMax, yMin, yMax, bXLabel):  
        numBins = 100
        
        n, bins, patches = axes.hist(beforeCorrection, numBins, density=1, facecolor='blue', alpha=0.25, label='before correction')
        n, bins, patches = axes.hist(afterCorrection, numBins, density=1, facecolor='green', alpha=0.25, label='after correction')
        
        y = mlab.normpdf(bins, np.mean(beforeCorrection), np.std(beforeCorrection))
        axes.plot(bins, y, 'b--')
        
        y = mlab.normpdf(bins, np.mean(afterCorrection), np.std(afterCorrection))
        axes.plot(bins, y, 'g--')
        
        axes.axvline(ls='--', c='red')
        
        axes.set_xlabel(bXLabel)
        axes.set_ylabel('% test genomes')
        axes.set_xlim([xMin, xMax])
        axes.set_ylim([yMin, yMax])
        
        legend = axes.legend(loc='best')
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

if __name__ == '__main__':
    histogram = Histogram()
    histogram.plot()
    histogram.savePlot('./simulations/correct_comp_cont_0.70-0.20.cont_per_by_g2.comp_by_is.png')
