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

class MarkerGeneFinder():
    """This class runs prodigal and hmmer creating data for parsing"""
    def __init__(self, threads):  
        self.logger = logging.getLogger()                  
        self.totalThreads = threads

    def find(self, binFiles, outDir, markerFile):
        """Identify marker genes in each bin using prodigal and HMMER."""
        makeSurePathExists(outDir)
      
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

        calcProc = [mp.Process(target = self.__processBin, args = (outDir, markerFile, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportProgress, args = (len(binFiles), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put(None)
        writeProc.join()
              
    def __processBin(self, outDir, markerFile, queueIn, queueOut):
        """Thread safe bin processing."""      
        while True:    
            binFile = queueIn.get(block=True, timeout=None) 
            if binFile == None:
                break   
            
            binId = binIdFromFilename(binFile)
            outDir = os.path.join(outDir, binId)
            makeSurePathExists(outDir)
    
            # run Prodigal      
            prodigal = ProdigalRunner() 
            aaGeneFile = prodigal.run(binFile, outDir)
    
            # run HMMER
            hmmer = HMMERRunner()
            hmmer.search(markerFile, aaGeneFile, outDir, '--cpu ' + str(self.threadsPerSearch) + ' --notextw -E 1 --domE 1')
    
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
            