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
import logging

import defaultValues

from hmmer import HMMERModelExtractor
from hmmerModelParser import HmmModelParser

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
    
    def numMarkers(self):
        return self.size()[0]
    
    def numSets(self):
        return len(self.markerSet)
    
    def getMarkerGenes(self):
        """Get marker genes within marker set."""
        markerGenes = set()
        for m in self.markerSet:
            for marker in m:
                markerGenes.add(marker)
                
        return markerGenes
    
    def genomeCheck(self, hits, bIndividualMarkers):
        """Calculate genome completeness and contamination."""
        if bIndividualMarkers:
            multiCopyCount = 0
            for marker in self.getMarkerGenes():
                if marker in hits:
                    multiCopyCount += (len(hits[marker]) - 1)
            
            percComp = 100 * float(len(hits)) / self.numMarkers()
            percCont = 100 * float(multiCopyCount) / self.numMarkers()
        else: 
            comp = 0.0
            cont = 0.0    
            for ms in self.markerSet:
                present = 0
                multiCopy = 0
                for marker in ms:
                    count = len(hits.get(marker, []))
                    if count == 1:
                        present += 1
                    elif count > 1:
                        present += 1
                        multiCopy += (count-1)
    
                comp += float(present) / len(ms)
                cont += float(multiCopy) / len(ms)

            percComp = 100 * comp / len(self.markerSet)
            percCont = 100 * cont / len(self.markerSet)
        
        return percComp, percCont
     
class MarkerSetParser():
    """Parse marker set file."""
    
    def __init__(self, threads=1):
        self.logger = logging.getLogger()
        self.numThreads = threads
        
        self.TAXONOMIC_MARKER_FILE = 1
        self.TREE_MARKER_FILE = 2
        self.HMM_MODELs_FILE = 3
        
    def getMarkerSets(self, outDir, binIds, markerFile):
        """Determine marker set for each bin."""
        
        # determine type of marker set file
        markerFileType = self.__markerFileType(markerFile)
        
        # get marker set for each bin
        binIdToMarkerSet = {}
        
        if markerFileType == self.TAXONOMIC_MARKER_FILE:
            markerSet = self.__parseTaxonomicMarkerSetFile(markerFile)
            
            for binId in binIds:
                binIdToMarkerSet[binId] = markerSet
        elif markerFileType == self.TREE_MARKER_FILE:
            binIdToMarkerSet = self.__parseLineageMarkerSetFile(markerFile)
        else:
            markerSet = [set()]
            modelParser = HmmModelParser(markerFile)
            for model in modelParser.parse():
                markerSet[0].add(model.acc)
                
            for binId in binIds:
                binIdToMarkerSet[binId] = MarkerSet(markerSet)
        
        return binIdToMarkerSet
            
    def createHmmModelFiles(self, outDir, binIds, markerFile):
        """Create HMM file for each bins marker set."""

        # determine type of marker set file
        markerFileType = self.__markerFileType(markerFile)

        # get HMM file for each bin
        binIdToHmmModelFile = {}
        if markerFileType == self.TAXONOMIC_MARKER_FILE:
            markerSet = self.__parseTaxonomicMarkerSetFile(markerFile)
            hmmModelFile = os.path.join(outDir, 'storage', 'hmms', 'taxonomic.hmm')
            self.__createMarkerHMMs(markerSet, hmmModelFile)
            
            for binId in binIds:
                binIdToHmmModelFile[binId] = hmmModelFile
        elif markerFileType == self.TREE_MARKER_FILE:
            binIdToMarkerSet = self.__parseLineageMarkerSetFile(markerFile)

            self.logger.info('  Extracting lineage-specific HMMs with %d threads:' % self.numThreads)
            for i, binId in enumerate(binIdToMarkerSet):
                if self.logger.getEffectiveLevel() <= logging.INFO:
                    statusStr = '    Finished extracting HMMs for %d of %d (%.2f%%) bins.' % (i+1, len(binIdToMarkerSet), float(i+1)*100/len(binIdToMarkerSet))
                    sys.stdout.write('%s\r' % statusStr)
                    sys.stdout.flush()
                
                hmmModelFile = os.path.join(outDir, 'storage', 'hmms', binId + '.hmm')
                self.__createMarkerHMMs(binIdToMarkerSet[binId], hmmModelFile, False)
                binIdToHmmModelFile[binId] = hmmModelFile 
                
            if self.logger.getEffectiveLevel() <= logging.INFO:
                sys.stdout.write('\n')
        else:
            for binId in binIds:
                binIdToHmmModelFile[binId] = markerFile
                
        return binIdToHmmModelFile
                
    def __markerFileType(self, markerFile):
        """Determine type of marker file."""
        with open(markerFile, 'r') as f:
            header = f.readline()
            
        if defaultValues.TAXON_MARKER_FILE_HEADER in header:
            return self.TAXONOMIC_MARKER_FILE
        elif defaultValues.LINEAGE_MARKER_FILE_HEADER in header:
            return self.TREE_MARKER_FILE
        elif 'HMMER3' in header:
            return self.HMM_MODELs_FILE
        else:
            self.logger.error('Unrecognized file type: ' + markerFile)
            sys.exit()
            
    def __createMarkerHMMs(self, markerSet, outputFile, bReportProgress = True):
        """Create HMM file for taxonomic markers."""

        # get list of marker genes  
        markersGene = markerSet.getMarkerGenes()
        
        # extract marker genes
        modelExtractor = HMMERModelExtractor(self.numThreads)
        modelExtractor.extract(defaultValues.HMM_MODELS, markersGene, outputFile, bReportProgress)
              
    def __parseTaxonomicMarkerSetFile(self, markerSetFile):
        """Parse marker set from a taxonomic-specific marker set file."""
        with open(markerSetFile) as f:
            f.readline() # skip header
            
            lineSplit = f.readline().split('\t')
            markerSet = MarkerSet(eval(lineSplit[4].rstrip()))
        
        return markerSet
    
    def __parseLineageMarkerSetFile(self, markerSetFile):
        """Parse marker set from a lineage-specific marker set file."""
        binIdToMarkerSets = {}
        with open(markerSetFile) as f:
            f.readline() # skip header
            
            for line in f:
                lineSplit = line.split('\t')
                binId = lineSplit[0]
                binIdToMarkerSets[binId] = MarkerSet(eval(lineSplit[5].rstrip()))
        
        return binIdToMarkerSets   