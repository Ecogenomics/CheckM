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

from lib.pfam import PFAM

class BinMarkerSets():
    """A collection of one or more marker sets associated with a bin."""
    def __init__(self, binId):
        self.logger = logging.getLogger()
        self.markerSets = []
        self.binId = binId
             
    def numMarkerSets(self):
        """Number of marker sets associated with bin."""
        return len(self.markerSets)
    
    def addMarkerSet(self, markerSet):
        """Add marker set to bin."""
        self.markerSets.append(markerSet)
        
    def markerSetIter(self):
        """Generator function for iterating over marker sets."""
        for markerSet in self.markerSets:
            yield markerSet      
            
    def getMarkerGenes(self):
        """Get marker genes from all marker sets."""
        markerGenes = set()
        for ms in self.markerSets:
            markerGenes.update(ms.getMarkerGenes())
                
        return markerGenes
            
    def write(self, fout):
        """Write marker set to file."""
        fout.write(self.binId)
        fout.write('\t' + str(len(self.markerSets)))
        for ms in self.markerSets:
            fout.write('\t' + str(ms)) 
        fout.write('\n')
        
    def read(self, line):
        """Construct bin marker set data from line."""
        lineSplit = line.split('\t')
        numMarkerSets = int(lineSplit[1])
        for i in xrange(0, numMarkerSets):
            uid = lineSplit[i*4+2]
            lineageStr = lineSplit[i*4+3]
            numGenomes = int(lineSplit[i*4+4])
            markerSet = eval(lineSplit[i*4+5])
            self.markerSets.append(MarkerSet(uid, lineageStr, numGenomes, markerSet))

class MarkerSet():
    """A collection of marker genes organized into co-located sets."""
    def __init__(self,  UID, lineageStr, numGenomes, markerSet):
        self.logger = logging.getLogger()
        
        self.UID = UID                  # unique ID of marker set
        self.lineageStr = lineageStr    # taxonomic string associated with marker set
        self.numGenomes = numGenomes    # number of genomes used to calculate marker set
        self.markerSet = markerSet      # marker genes organized into co-located sets
           
    def __repr__(self):
        return str(self.UID) + '\t' + self.lineageStr + '\t' + str(self.numGenomes) + '\t' + str(self.markerSet)
        
    def size(self):
        """Number of marker genes and marker gene sets."""
        numMarkerGenes = 0
        for m in self.markerSet:
            numMarkerGenes += len(m)
            
        return numMarkerGenes, len(self.markerSet)
        
    def numMarkers(self):
        """Number of marker genes."""
        return self.size()[0]

    def numSets(self):
        """Number of marker sets."""
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
            present = 0
            multiCopyCount = 0
            for marker in self.getMarkerGenes():
                if marker in hits:
                    present += 1
                    multiCopyCount += (len(hits[marker]) - 1)
            
            percComp = 100 * float(present) / self.numMarkers()
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
        binIdToBinMarkerSets = {}
        
        if markerFileType == self.TAXONOMIC_MARKER_FILE:
            binMarkerSets = self.__parseTaxonomicMarkerSetFile(markerFile)
            
            for binId in binIds:
                binIdToBinMarkerSets[binId] = binMarkerSets
        elif markerFileType == self.TREE_MARKER_FILE:
            binIdToBinMarkerSets = self.__parseLineageMarkerSetFile(markerFile)
        else:
            markers = [set()]
            modelParser = HmmModelParser(markerFile)
            for model in modelParser.parse():
                markers[0].add(model.acc)
            markerSet = MarkerSet(0, "N/A", -1, markers)
            
            for binId in binIds:
                binMarkerSets = BinMarkerSets(binId)
                binMarkerSets.addMarkerSet(markerSet)
                binIdToBinMarkerSets[binId] = binMarkerSets
        
        return binIdToBinMarkerSets
            
    def createHmmModelFiles(self, outDir, binIds, markerFile):
        """Create HMM file for each bins marker set."""

        # determine type of marker set file
        markerFileType = self.__markerFileType(markerFile)

        # get HMM file for each bin
        binIdToHmmModelFile = {}
        if markerFileType == self.TAXONOMIC_MARKER_FILE:
            binMarkerSets = self.__parseTaxonomicMarkerSetFile(markerFile)
            hmmModelFile = os.path.join(outDir, 'storage', 'hmms', 'taxonomic.hmm')
            self.__createMarkerHMMs(binMarkerSets, hmmModelFile)
            
            for binId in binIds:
                binIdToHmmModelFile[binId] = hmmModelFile
        elif markerFileType == self.TREE_MARKER_FILE:
            binIdToBinMarkerSets = self.__parseLineageMarkerSetFile(markerFile)

            self.logger.info('  Extracting lineage-specific HMMs with %d threads:' % self.numThreads)
            for i, binId in enumerate(binIdToBinMarkerSets):
                if self.logger.getEffectiveLevel() <= logging.INFO:
                    statusStr = '    Finished extracting HMMs for %d of %d (%.2f%%) bins.' % (i+1, len(binIdToBinMarkerSets), float(i+1)*100/len(binIdToBinMarkerSets))
                    sys.stdout.write('%s\r' % statusStr)
                    sys.stdout.flush()
                
                hmmModelFile = os.path.join(outDir, 'storage', 'hmms', binId + '.hmm')
                self.__createMarkerHMMs(binIdToBinMarkerSets[binId], hmmModelFile, False)
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
            
    def __createMarkerHMMs(self, binMarkerSet, outputFile, bReportProgress = True):
        """Create HMM file for taxonomic markers."""

        # get list of marker genes  
        markerGenes = binMarkerSet.getMarkerGenes()
        
        # get all genes from the same clan as any marker gene
        pfam = PFAM()
        genesInSameClan = pfam.genesInSameClan(markerGenes)
        
        # extract marker genes along with all genes from the same clan
        allMarkers = markerGenes | genesInSameClan
        
        if bReportProgress:
            self.logger.info("  There are %d genes in the marker set and %d genes from the same PFAM clan." % (len(markerGenes), len(genesInSameClan)))
            
        modelExtractor = HMMERModelExtractor(self.numThreads)
        modelExtractor.extract(defaultValues.HMM_MODELS, allMarkers, outputFile, bReportProgress)
              
    def __parseTaxonomicMarkerSetFile(self, markerSetFile):
        """Parse marker set from a taxonomic-specific marker set file."""
        with open(markerSetFile) as f:
            f.readline() # skip header
            
            binLine = f.readline()
            taxonId = binLine.split('\t')[0]
            binMarkerSets = BinMarkerSets(taxonId)
            binMarkerSets.read(binLine)
        
        return binMarkerSets
    
    def __parseLineageMarkerSetFile(self, markerSetFile):
        """Parse marker sets from a lineage-specific marker set file."""
        binIdToBinMarkerSets = {}
        with open(markerSetFile) as f:
            f.readline() # skip header
             
            for line in f:
                lineSplit = line.split('\t')
                binId = lineSplit[0]
                
                binMarkerSets = BinMarkerSets(binId)
                binMarkerSets.read(line)
                
                binIdToBinMarkerSets[binId] = binMarkerSets
        
        return binIdToBinMarkerSets   