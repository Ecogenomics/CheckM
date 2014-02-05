###############################################################################
#
# markerGeneFinder.py - identify marker genes in genome bins
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
import multiprocessing as mp
import logging

from common import binIdFromFilename, makeSurePathExists

import defaultValues

from hmmer import HMMERRunner, HMMERModelExtractor
from prodigal import ProdigalRunner

from markerSet import parseTaxonomicMarkerSetFile, parseLineageMarkerSetFile

class MarkerGeneFinder():
    """Identify marker genes within binned sequences using Prodigal and HMMER."""
    def __init__(self, threads):  
        self.logger = logging.getLogger()                  
        self.totalThreads = threads
        
        self.TAXONOMIC_MARKER_FILE = 1
        self.TREE_MARKER_FILE = 2
        self.HMM_MODELs_FILE = 3

    def find(self, binFiles, outDir, tableOut, hmmerOut, markerFile):
        """Identify marker genes in each bin using prodigal and HMMER."""
        
        # determine type of marker file
        markerFileType = self.__markerFileType(markerFile)
        
        # get list of all bins
        binIds = []
        for binFile in binFiles:
            binIds.append(binIdFromFilename(binFile))
        
        # create HMM model file(s)
        makeSurePathExists(os.path.join(outDir, 'storage', 'tmp'))
        
        hmmModelFiles = {}
        if markerFileType == self.TAXONOMIC_MARKER_FILE:
            markerSet = parseTaxonomicMarkerSetFile(markerFile)
            hmmModelFile = os.path.join(outDir, 'storage', 'tmp', 'taxonomic.hmm')
            self.__createMarkerHMMs(markerSet, hmmModelFile)
            
            for binId in binIds:
                hmmModelFiles[binId] = hmmModelFile
        elif markerFileType == self.TREE_MARKER_FILE:
            binIdToMarkerSets = parseLineageMarkerSetFile(markerFile)

            self.logger.info('  Extracting lineage-specific HMM models with %d threads:' % self.totalThreads)
            for i, binId in enumerate(binIdToMarkerSets):
                if self.logger.getEffectiveLevel() <= logging.INFO:
                    statusStr = '    Finished extracting models from %d of %d (%.2f%%) bins.' % (i, len(binIdToMarkerSets), float(i)*100/len(binIdToMarkerSets))
                    sys.stdout.write('%s\r' % statusStr)
                    sys.stdout.flush()
                
                hmmModelFile = os.path.join(outDir, 'storage', 'tmp', binId + '.hmm')
                self.__createMarkerHMMs(binIdToMarkerSets[binId], hmmModelFile, False)
                hmmModelFiles[binId] = hmmModelFile 
                
            if self.logger.getEffectiveLevel() <= logging.INFO:
                sys.stdout.write('\n')
        else:
            for binId in binIds:
                hmmModelFiles[binId] = markerFile

        # process each fasta file
        self.threadsPerSearch = max(1, int(self.totalThreads / len(binFiles)))
        self.logger.info("  Identifying marker genes in %d bins with %d threads:" % (len(binFiles), self.totalThreads))
        
        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for binFile in binFiles:
            workerQueue.put(binFile)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__processBin, args = (outDir, tableOut, hmmerOut, hmmModelFiles, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportProgress, args = (len(binFiles), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put(None)
        writeProc.join()
        
        return hmmModelFiles
        
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
        modelExtractor = HMMERModelExtractor(self.totalThreads)
        modelExtractor.extract(defaultValues.HMM_MODELS, markersGene, outputFile, bReportProgress)
              
    def __processBin(self, outDir, tableOut, hmmerOut, hmmModelFiles, queueIn, queueOut):
        """Thread safe bin processing."""      
        while True:    
            binFile = queueIn.get(block=True, timeout=None) 
            if binFile == None:
                break   
            
            binId = binIdFromFilename(binFile)
            binDir = os.path.join(outDir, binId)
            makeSurePathExists(binDir)
    
            # run Prodigal     
            prodigal = ProdigalRunner(binDir) 
            if not prodigal.areORFsCalled():
                prodigal.run(binFile)
    
            # run HMMER
            hmmer = HMMERRunner()
            tableOutPath = os.path.join(binDir, tableOut)
            hmmerOutPath = os.path.join(binDir, hmmerOut)
            hmmer.search(hmmModelFiles[binId], prodigal.aaGeneFile, tableOutPath, hmmerOutPath, '--cpu ' + str(self.threadsPerSearch) + ' --notextw -E 1 --domE 1')
    
            queueOut.put(binId)
            
    def __reportProgress(self, numBins, queueIn):
        """Report number of processed bins."""      
        
        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        
        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
         
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
            