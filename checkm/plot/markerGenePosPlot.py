###############################################################################
#
# markerGenePosPlot.py - Create a plot showing the position of marker genes.
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

import os
import sys

from AbstractPlot import AbstractPlot

from checkm.seqUtils import readFasta

import checkm.defaultValues as defaultValues
from checkm.hmmer import HMMERParser

class markerGenePosPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)
        
    def getMarkerGenePositions(self, fastaFile, resultsDir):
        markerGenePositions = {}
        
        inputFile = resultsDir + '/' + fastaFile + '/' + defaultValues.__CHECKM_DEFAULT_HMMER_TXT_OUT__
        if not os.path.exists(inputFile):
            sys.stderr.write('[Error] Expected input file does not exists: ' + inputFile)
            sys.exit()
        
        with open(inputFile, 'r') as hmmer_handle:
                try:
                    HP = HMMERParser(hmmer_handle)
                except:
                    sys.stderr.write("Error opening HMM file: ", inputFile)
                    raise
                
                while True:
                    hit = HP.next()
                    if hit is None:
                        break
                    
                    seqId = hit.target_name[0:hit.target_name.rfind('_')]
                    markerGenePositions[seqId] = markerGenePositions.get(seqId, []) + [hit.ali_to]
        
        return markerGenePositions
    
    def plot(self, fastaFile, resultsDir):
        # Get length of each sequence in bin
        seqs = readFasta(fastaFile)
        seqLens = {}
        for seqId, seq in seqs.iteritems():
            seqLens[seqId] = len(seq)
            
        # Get position of marker genes
        markerGenePositions = self.getMarkerGenePositions(fastaFile, resultsDir)
        
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        heightBottomLabels = 0.3 + self.options.fig_padding   # inches
        widthSideLabel = 0.5 + self.options.fig_padding       # inches 
        axes = self.fig.add_axes([widthSideLabel/self.options.width,heightBottomLabels/self.options.height,\
                                                                        1.0-(widthSideLabel+self.options.fig_padding )/self.options.width,\
                                                                        1.0-(heightBottomLabels+self.options.fig_padding )/self.options.height])
        
            
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
                          
        self.draw()