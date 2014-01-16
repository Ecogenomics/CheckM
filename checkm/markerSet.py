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

import os
import sys

from lib.img import IMG

from common import makeSurePathExists

class MarkerSet():
    """Calculate and process marker sets."""
    def __init__(self):
        pass
    
    def identifyLineageSpecificMarkers(self, treeFile, binFile):
        """Calculate lineage-specific marker set for a given genome bin."""
        
        # determine lineage suitable for calculating markers
        genomeIds = []
          
        # build gene count table
        img = IMG()
        countTable = img.countTable(genomeIds)
        countTable = img.filterTable(genomeIds, countTable, 0.9*ubiquityThreshold, 0.9*singleCopyThreshold)

        # identify marker genes for genomes
        markerGenes = markerset.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
        tigrToRemove = img.identifyRedundantTIGRFAMs(markerGenes)
        markerGenes = markerGenes - tigrToRemove

        # identify marker sets
        geneDistTable = img.geneDistTable(genomeIds, markerGenes)
        colocatedGenes = markerset.colocatedGenes(geneDistTable)
        colocatedSets = markerset.colocatedSets(colocatedGenes, markerGenes)

   
    def checkBin(self):
        pass