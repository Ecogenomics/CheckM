###############################################################################
#
# coverage.py - calculate coverage of all sequences
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

import sys
import os
import multiprocessing as mp
import collections

import pysam

import defaultValues

from common import reassignStdOut, restoreStdOut, binIdFromFilename
from seqUtils import readFasta

class CoverageStruct():
    def __init__(self, seqLen, mappedReads, coverage):
        self.seqLen = seqLen
        self.mappedReads = mappedReads
        self.coverage = coverage

class Coverage():
    """Calculate coverage of all sequences."""
    def __init__(self, threads=1):
        # thready stuff
        self.totalThreads = threads
        self.varLock = mp.Lock()
        self.threadPool = mp.BoundedSemaphore(self.totalThreads)
        
        self.numFilesStarted = mp.Value('i', 0)
        self.numFiles = 0
        
        self.CoverageStruct = collections.namedtuple("CoverageInfo", "seqLen mappedReads coverage")

    def calculate(self, binFiles, bamFiles, outFile, bPairsOnly=False, bQuiet=False):
        """Calculate coverage of sequences for each BAM file."""
        
        # determine bin assignment of each sequence
        if not bQuiet:
            print '  Determining bin assignment of each sequence.'
            
        seqIdToBinId = {}
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)
            
            seqs = readFasta(binFile)
            for seqId in seqs.keys():
                seqIdToBinId[seqId] = binId
        
        # process each fasta file
        if not bQuiet:
            print "  Processing %d file(s) with %d threads:" % (len(bamFiles), self.totalThreads)
            
        # make sure all BAM files are sorted
        self.numFiles = len(bamFiles)
        for bamFile in bamFiles: 
            if not os.path.exists(bamFile + '.bai'):
                sys.stderr.write('[Error] BAM file is not sorted: ' + bamFile + '\n')
                sys.exit()
 
        # calculate coverage of each BAM file
        processes = []
        coverageInfo = {}
        for bamFile in bamFiles: 
            coverageInfo[bamFile] = mp.Manager().dict()
            p = mp.Process(target=self.__processBam, args=(bamFile, coverageInfo[bamFile], bPairsOnly, bQuiet))
            p.start()
            processes.append(p)
            
        for p in processes:
            p.join()
            
        # redirect output
        oldStdOut = reassignStdOut(outFile)
        
        header = 'Sequence Id\tBin Id\tSequence length (bp)\t'
        for bamFile in bamFiles:
            header += 'Bam Id\tCoverage\tMapped reads'
        
        print(header)

        for seqId in coverageInfo[coverageInfo.keys()[0]].keys():
            rowStr = seqId + '\t' + seqIdToBinId.get(seqId, defaultValues.__CHECKM_DEFAULT_UNBINNED__) + '\t' + str(coverageInfo[coverageInfo.keys()[0]][seqId].seqLen)
            for bamFile in bamFiles:
                bamId = binIdFromFilename(bamFile)
                rowStr += '\t' + bamId + '\t' + str(coverageInfo[bamFile][seqId].coverage) + '\t' + str(coverageInfo[bamFile][seqId].mappedReads)
            print(rowStr)
        
        # restore stdout
        restoreStdOut(outFile, oldStdOut, bQuiet)  

    def __processBam(self, bamFile, coverageInfo, bPairsOnly, bQuiet):
        """Calculate coverage of sequences in BAM file."""
        self.threadPool.acquire()
        try:
            self.varLock.acquire()
            try:
                self.numFilesStarted.value += 1
                if not bQuiet:
                    print '    Processing %s (%d of %d)' % (bamFile, self.numFilesStarted.value, self.numFiles)
            finally:
                self.varLock.release()

            bamIn = pysam.Samfile(bamFile, 'rb')
            
            refLen = {}
            for refName, length in zip(bamIn.references, bamIn.lengths):
                refLen[refName] = length
    
            # read entire sorted BAM file and determine coverage for each reference sequence
            mappedReads = 0
            coverage = 0
            curRef = None
            for read in bamIn.fetch(until_eof = True):
                if read.is_unmapped or read.is_duplicate or read.is_qcfail or read.is_secondary:
                    continue
    
                if bPairsOnly and not read.is_proper_pair:
                    continue
                
                if curRef == None:
                    curRef = bamIn.getrname(read.tid)
                    
                refName = bamIn.getrname(read.tid)
                if curRef == refName:
                    mappedReads += 1
                    coverage += read.rlen
                else:
                    refLength = refLen[curRef]
                    coverage = float(coverage) / refLength
                    coverageInfo[curRef] = CoverageStruct(seqLen = refLength, mappedReads = mappedReads, coverage = coverage)
                    
                    mappedReads = 1
                    coverage = read.rlen
                    curRef = refName
                
            # results for final reference sequence    
            refLength = refLen[curRef]
            coverage /= refLength
            coverageInfo[curRef] = CoverageStruct(seqLen = refLength, mappedReads = mappedReads, coverage = coverage)

            bamIn.close()
            
        finally:
            self.threadPool.release()

    def parseCoverage(self, coverageFile):
        """Read coverage information from file."""
        coverageStats = {}
        bHeader = True
        for line in open(coverageFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            binId = lineSplit[1]
            #seqLen = lineSplit[2]
            
            if binId not in coverageStats:
                coverageStats[binId] = {}
                    
            if seqId not in coverageStats[binId]:
                coverageStats[binId][seqId] = {}
                
            for i in xrange(3, len(lineSplit), 3):
                bamId = lineSplit[i]
                coverage = float(lineSplit[i+1])
                #mappedReads = int(lineSplit[i+2])
                coverageStats[binId][seqId][bamId] = coverage
                
        return coverageStats