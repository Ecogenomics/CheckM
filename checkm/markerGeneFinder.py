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

from os import mkdir
import os.path
import threading
import time

from common import makeSurePathExists

from hmmer import HMMERRunner
from prodigal import ProdigalRunner

class MarkerGeneFinder():
    """This class runs prodigal and hmmer creating data for parsing"""
    def __init__(self, threads=1):         
        # thready stuff            
        self.varLock = threading.Lock() # we don't have many variables here, so use one lock for everything    
        self.totalThreads = threads
        self.threadsPerSearch = 1
        self.threadPool = threading.BoundedSemaphore(self.totalThreads)

        self.numFiles = 0
        self.numFilesStarted = 0
        self.numFilesParsed = 0

    def find(self, inFiles, outFolder, hmm, bQuiet=False):
        """Identify marker genes in each input fasta file using prodigal and HMMER."""
        if not os.path.exists(outFolder):
            mkdir(outFolder)
            
        self.hmmer = HMMERRunner()
        self.prodigal = ProdigalRunner()
            
        # process each fasta file
        self.numFiles = len(inFiles)
        self.threadsPerSearch = max(1, int(self.totalThreads / self.numFiles))
        if not bQuiet:
            print "  Processing %d files with %d threads:" % (self.numFiles, self.totalThreads)
        
        for fasta in inFiles:
            t = threading.Thread(target=self.__processFasta, args=(fasta, outFolder, hmm, bQuiet))
            t.start()
            
        while True:
            # don't exit till we're done!
            self.varLock.acquire()
            doned = False
            try:
                doned = self.numFilesParsed >= self.numFiles
            finally:
                self.varLock.release()
            
            if doned:
                break
            else:
                time.sleep(1)
              
    def __processFasta(self, fasta, outFolder, hmm, bQuiet=False):
        """Thread safe fasta processing"""
        self.threadPool.acquire()
        try:
            self.varLock.acquire()
            try:
                # create output directory for fasta file
                outDir = os.path.join(outFolder, os.path.basename(fasta))
                makeSurePathExists(outDir)
                
                self.numFilesStarted += 1
                if not bQuiet:
                    print "    Processing %s (%d of %d)" % (fasta, self.numFilesStarted, self.numFiles) 
            finally:
                self.varLock.release()

            # run Prodigal       
            aaGeneFile = self.prodigal.run(fasta, outDir)

            # run HMMER
            self.hmmer.search(hmm, aaGeneFile, outDir, '--cpu ' + str(self.threadsPerSearch))

            # let the world know we've parsed this file
            self.varLock.acquire()
            try:
                self.numFilesParsed += 1
            finally:
                self.varLock.release()
            
        finally:
            self.threadPool.release()
