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

from hmmer import HMMERRunner
from prodigal import ProdigalRunner

from markerSets import MarkerSetParser

class MarkerGeneFinder():
    """Identify marker genes within binned sequences using Prodigal and HMMER."""
    def __init__(self, threads):  
        self.logger = logging.getLogger()                  
        self.totalThreads = threads

    def find(self, binFiles, outDir, tableOut, hmmerOut, markerFile, bKeepAlignment, bNucORFs):
        """Identify marker genes in each bin using prodigal and HMMER."""
                
        # get bin ids
        binIds = []
        for binFile in binFiles:
            binIds.append(binIdFromFilename(binFile))
                
        # create HMM file(s)
        markerSetParser = MarkerSetParser(self.totalThreads)
        binIdToHmmModelFile = markerSetParser.createHmmModelFiles(outDir, binIds, markerFile)

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

        calcProc = [mp.Process(target = self.__processBin, args = (outDir, tableOut, hmmerOut, binIdToHmmModelFile, bKeepAlignment, bNucORFs, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportProgress, args = (len(binFiles), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put(None)
        writeProc.join()
        
        return binIdToHmmModelFile
              
    def __processBin(self, outDir, tableOut, hmmerOut, binIdToHmmModelFile, bKeepAlignment, bNucORFs, queueIn, queueOut):
        """Thread safe bin processing."""      
        while True:    
            binFile = queueIn.get(block=True, timeout=None) 
            if binFile == None:
                break   
            
            binId = binIdFromFilename(binFile)
            binDir = os.path.join(outDir, 'bins', binId)
            makeSurePathExists(binDir)
    
            # run Prodigal     
            prodigal = ProdigalRunner(binDir) 
            if not prodigal.areORFsCalled():
                prodigal.run(binFile, bNucORFs)
    
            # run HMMER
            hmmer = HMMERRunner()
            tableOutPath = os.path.join(binDir, tableOut)
            hmmerOutPath = os.path.join(binDir, hmmerOut)
            
            keepAlignStr = ''
            if not bKeepAlignment:
                keepAlignStr = '--noali'
            hmmer.search(binIdToHmmModelFile[binId], prodigal.aaGeneFile, tableOutPath, hmmerOutPath, '--cpu ' + str(self.threadsPerSearch) + ' --notextw -E 0.1 --domE 0.1 ' + keepAlignStr, bKeepAlignment)
    
            queueOut.put(binId)
            
    def __reportProgress(self, numBins, queueIn):
        """Report number of processed bins."""      
        
        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()
        
        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()
         
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')
            