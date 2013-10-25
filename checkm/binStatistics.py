###############################################################################
#
# binStatistics.py - calculate statistics for each putative genome bin
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
import threading
import time
import math

import defaultValues
from common import readFasta
from mathHelper import mean

class BinStatistics():
    """Calculate statistics (GC, coding density, etc.) for putative genome bins."""
    def __init__(self, threads=1):         
        # thready stuff            
        self.varLock = threading.Lock() # we don't have many variables here, so use one lock for everything    
        self.totalThreads = threads
        self.threadsPerSearch = 1
        self.threadPool = threading.BoundedSemaphore(self.totalThreads)

        self.numFiles = 0
        self.numFilesStarted = 0
        self.numFilesParsed = 0

    def calculate(self, inFiles, outFolder, quiet=False):
        """Calculate statistics for each putative genome bin."""

        # process each fasta file
        self.numFiles = len(inFiles)
        if not quiet:
            print "  Processing %d files with %d threads:" % (self.numFiles, self.totalThreads)
        
        binStats = {}
        scaffoldStats = {}
        for fasta in inFiles:
            binName = os.path.basename(fasta)
            binStats[binName] = {}
            scaffoldStats[binName] = {}
            
            t = threading.Thread(target=self.__processBin, args=(fasta, outFolder, binStats[binName], scaffoldStats[binName], quiet))
            t.start()
            
        while True:
            self.varLock.acquire()
            bDone = False
            try:
                bDone = self.numFilesParsed >= self.numFiles
            finally:
                self.varLock.release()
            
            if bDone:
                break
            else:
                time.sleep(1)
                
        fout = open(os.path.join(outFolder, defaultValues.__CHECKM_DEFAULT_BIN_STATS_FILE__), 'w')
        fout.write(str(binStats))
        fout.close()
        
        fout = open(os.path.join(outFolder, defaultValues.__CHECKM_DEFAULT_SCAFFOLD_STATS_FILE__), 'w')
        fout.write(str(scaffoldStats))
        fout.close()
              
    def __processBin(self, fasta, outFolder, binStats, scaffoldStats, quiet=False):
        """Thread safe fasta processing"""
        self.threadPool.acquire()
        try:
            self.varLock.acquire()
            try:
                self.numFilesStarted += 1
                if not quiet:
                    print '    Processing %s (%d of %d)' % (fasta, self.numFilesStarted, self.numFiles) 
            finally:
                self.varLock.release()

            outDir = os.path.join(outFolder, os.path.basename(fasta))

            # read seqs
            seqs = readFasta(fasta)
            for seqId in seqs:
                scaffoldStats[seqId] = {}
                
            # calculate GC statistics
            GC, stdGC = self.calculateGC(seqs, scaffoldStats) 
            binStats['GC'] = GC
            binStats['GC std'] = stdGC
            
            # calculate statistics related to scaffold lengths
            minSeqLen, maxSeqLen, genomeSize, N50, numContigs  = self.calculateScaffoldLengthStats(seqs, scaffoldStats)  
            binStats['# scaffolds'] = len(seqs)
            binStats['# contigs'] = numContigs
            binStats['Shortest scaffold'] = minSeqLen
            binStats['Longest scaffold'] = maxSeqLen
            binStats['Genome size'] = genomeSize
            binStats['N50'] = N50
            
            # calculate coding density statistics
            codingDensity, numORFs = self.calculateCodingDensity(outDir, genomeSize)
            binStats['Coding density'] = codingDensity
            binStats['# predicted ORFs'] = numORFs

            # let the world know we've parsed this file
            self.varLock.acquire()
            try:
                self.numFilesParsed += 1
            finally:
                self.varLock.release()
            
        finally:
            self.threadPool.release()
            
    def calculateGC(self, seqs, scaffoldStats):
        """Calculate fraction of nucleotides that are G or C."""
        totalGC = 0
        totalAT = 0
        gcPerSeq = []
        for seqId, seq in seqs.iteritems():  
            testSeq = seq.upper()
            gc = testSeq.count('G') + testSeq.count('C')
            at = testSeq.count('A') + testSeq.count('T') + testSeq.count('U')
            
            totalGC += gc
            totalAT += at
            
            gcContent = float(gc) / (gc + at) 
            scaffoldStats[seqId]['GC'] = gcContent
            
            if len(seq) > defaultValues.__CHECKM_DEFAULT_MIN_SEQ_LEN_GC_STD__:
                gcPerSeq.append(gcContent)
                
        GC = float(totalGC) / (totalGC + totalAT)
        
        varGC = mean(map(lambda x: (x - GC)**2, gcPerSeq))

        return GC, math.sqrt(varGC)
    
    def calculateScaffoldLengthStats(self, seqs, scaffoldStats):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        seqLens = []
        totalLen = 0
        numContigs = 0
        for seqId, seq in seqs.iteritems():
            seqLen = len(seq)
            seqLens.append(seqLen)
            totalLen += seqLen
            
            scaffoldStats[seqId]['length'] = seqLen
            
            splitSeq = seq.split(defaultValues.__CHECKM_DEFAULT_CONTIG_BREAK__)
            numContigs += len([x for x in splitSeq if len(x) > len(defaultValues.__CHECKM_DEFAULT_CONTIG_BREAK__)])
            
        thresholdN50 = totalLen / 2.0
            
        seqLens.sort(reverse=True)
        
        testSum = 0
        for seqLen in seqLens:
            testSum += seqLen
            if testSum >= thresholdN50:
                N50 = seqLen
                break
            
        return min(seqLens), max(seqLens), sum(seqLens), N50, numContigs
    
    def calculateCodingDensity(self, outDir, genomeSize):
        """Calculate coding density of putative genome bin."""
        ntFile = os.path.join(outDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_NT__)
        ntGenes = readFasta(ntFile)
        
        codingBasePairs = 0
        for _, gene in ntGenes.iteritems():
            codingBasePairs += len(gene)
            
        return float(codingBasePairs) / genomeSize, len(ntGenes)
        
            