###############################################################################
#
# gcPlots.py - Create a GC histogram and delta-GC plot. 
#
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

import matplotlib as mpl
import matplotlib.pyplot as pylab

import numpy as np

from AbstractPlot import AbstractPlot

from checkm.seqUtils import readFasta, baseCount

class gcPlots(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)
    
    def plot(self, fastaFile):
        # Global plot settings
        mpl.rcParams['font.size'] = self.options.font_size
        mpl.rcParams['axes.titlesize'] = self.options.font_size
        mpl.rcParams['axes.labelsize'] = self.options.font_size
        mpl.rcParams['xtick.labelsize'] = self.options.font_size
        mpl.rcParams['ytick.labelsize'] = self.options.font_size
        mpl.rcParams['legend.fontsize'] = self.options.font_size
        
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        heightBottomLabels = 0.3 + self.options.fig_padding   # inches
        widthSideLabel = 0.3 + self.options.fig_padding       # inches 
        axesHist = self.fig.add_axes([widthSideLabel/self.options.width,\
                                        heightBottomLabels/self.options.height,\
                                        0.5 - (widthSideLabel + self.options.fig_padding)/self.options.width,\
                                        1.0 - (heightBottomLabels + self.options.fig_padding)/self.options.height])
        
        axesDeltaGC = self.fig.add_axes([0.5 + widthSideLabel/self.options.width,\
                                            heightBottomLabels/self.options.height,\
                                            0.5 - (widthSideLabel + self.options.fig_padding)/self.options.width,\
                                            1.0 - (heightBottomLabels + self.options.fig_padding)/self.options.height])
        
        # get GC for windows
        seqs = readFasta(fastaFile)

        data = []
        for _, seq in seqs.iteritems():
            start = 0
            end = self.options.window_size
            while(end < len(seq)):
                a, c, g, t = baseCount(seq[start:end])
                data.append(float(g + c) / (a + c + g + t))
                
                start = end
                end += self.options.window_size
            
        # Histogram plot 
        bins = [0.0]
        binWidth = self.options.gc_bin_width
        binEnd = binWidth
        while binEnd <= 1.0:
            bins.append(binEnd)
            binEnd += binWidth
            
        axesHist.hist(data, bins=bins, normed=True, color=(0.5,0.5,0.5))    
        axesHist.set_xlabel('% GC')
        axesHist.set_ylabel('% windows (' + str(self.options.window_size) + ' bp)')
            
        # Prettify plot     
        for a in axesHist.yaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
                
        for a in axesHist.xaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
            
        for line in axesHist.yaxis.get_ticklines(): 
            line.set_color(self.axesColour)
                
        for line in axesHist.xaxis.get_ticklines(): 
            line.set_color(self.axesColour)
            
        for loc, spine in axesHist.spines.iteritems():
            if loc in ['right','top']:
                spine.set_color('none') 
            else:
                spine.set_color(self.axesColour)
                
        # get delta-GC and sequence lengths
        GCs = []
        seqLens = []
        gcTotal = 0
        basesTotal = 0
        for _, seq in seqs.iteritems():
            a, c, g, t = baseCount(seq)
            gc = g + c
            bases = a + c + g + t
            
            GCs.append(float(gc) / (bases))
            
            gcTotal += gc
            basesTotal += bases
            
            seqLens.append(len(seq))

        meanGC = float(gcTotal) / basesTotal
        deltaGCs = np.array(GCs) - meanGC
                
        # Delta-GC vs Sequence length plot 
        axesDeltaGC.scatter(deltaGCs, seqLens, c=abs(deltaGCs), s=10, lw=0.5, cmap=pylab.cm.Greys)    
        axesDeltaGC.set_xlabel(r'$\Delta$ GC')
        axesDeltaGC.set_ylabel('Sequence length (Kbps)')
        
        # ensure y-axis include zero
        _, yMax = axesDeltaGC.get_ylim()
        axesDeltaGC.set_ylim([0, yMax])
        
        # draw vertical line at x=0
        axesDeltaGC.vlines(0, 0, yMax, linestyle='dashed', color=self.axesColour, zorder=0)
        
        # Change sequence lengths from bps to kbps
        yticks = axesDeltaGC.get_yticks()
        kbpLabels = []
        for seqLen in yticks:
            label = '%.1f' % (float(seqLen)/1000)
            label = label.replace('.0', '') # remove trailing zero
            kbpLabels.append(label)
        axesDeltaGC.set_yticklabels(kbpLabels)
            
        # Prettify plot     
        for a in axesDeltaGC.yaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
                
        for a in axesDeltaGC.xaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
            
        for line in axesDeltaGC.yaxis.get_ticklines(): 
            line.set_color(self.axesColour)
                
        for line in axesDeltaGC.xaxis.get_ticklines(): 
            line.set_color(self.axesColour)
            
        for loc, spine in axesDeltaGC.spines.iteritems():
            if loc in ['right','top']:
                spine.set_color('none') 
            else:
                spine.set_color(self.axesColour)
                          
        self.draw()