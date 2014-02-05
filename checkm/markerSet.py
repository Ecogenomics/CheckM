###############################################################################
#
# markerSet.py - Calculate and process marker sets.
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

import logging

def parseTaxonomicMarkerSetFile(markerSetFile):
    """Parse marker set from a taxonomic-specific marker set file."""
    with open(markerSetFile) as f:
        f.readline() # skip header
        
        lineSplit = f.readline().split('\t')
        markerSet = MarkerSet(eval(lineSplit[4].rstrip()))
    
    return markerSet

def parseLineageMarkerSetFile(markerSetFile):
    """Parse marker set from a lineage-specific marker set file."""
    binIdToMarkerSets = {}
    with open(markerSetFile) as f:
        f.readline() # skip header
        
        for line in f:
            lineSplit = line.split('\t')
            binId = lineSplit[0]
            binIdToMarkerSets[binId] = MarkerSet(eval(lineSplit[4].rstrip()))
    
    return binIdToMarkerSets

class MarkerSet():
    """Marker genes organized into co-located sets."""
    def __init__(self, markerSet):
        """
          markerSet: a list of marker gene sets
        """
        self.logger = logging.getLogger()
        self.markerSet = markerSet
        
    def __repr__(self):
        return str(self.markerSet)
        
    def size(self):
        """Number of marker genes and marker gene sets."""
        numMarkerGenes = 0
        for m in self.markerSet:
            numMarkerGenes += len(m)
            
        return numMarkerGenes, len(self.markerSet)
    
    def getMarkerGenes(self):
        """Get marker genes within marker set."""
        markerGenes = set()
        for m in self.markerSet:
            for marker in m:
                markerGenes.add(marker)
                
        return markerGenes
        
    def identifyTaxonomicSpecificMarkers(self, rank, clade, 
                                            ubiquityThreshold, singleCopyThreshold,
                                            outputFile):
        """Calculate a taxonomic-specific marker sets."""
        pass

    def identifyLineageSpecificMarkers(self, treeFile, binFile):
        """Calculate lineage-specific marker set."""
        pass

    def checkBin(self):
        pass