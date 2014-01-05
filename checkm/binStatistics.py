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
import threading
import time
import math
import logging

import defaultValues
from seqUtils import readFasta, baseCount, calculateN50
from common import binIdFromFilename
from mathHelper import mean

class BinStatistics():
    """Calculate statistics (GC, coding density, etc.) for genome bins."""
    def __init__(self, threads):   
        self.logger = logging.getLogger()
              
        # thready stuff            
        self.varLock = threading.Lock() # we don't have many variables here, so use one lock for everything    
        self.totalThreads = threads
        self.threadsPerSearch = 1
        self.threadPool = threading.BoundedSemaphore(self.totalThreads)

        self.numFiles = 0
        self.numFilesStarted = 0
        self.numFilesParsed = 0

    def calculate(self, inFiles, outFolder):
        """Calculate statistics for each putative genome bin."""
        logger = logging.getLogger()
        
        # process each fasta file
        self.numFiles = len(inFiles)
        
        logger.info("  Processing %d files with %d threads:" % (self.numFiles, self.totalThreads))
        
        binStats = {}
        seqStats = {}
        for fasta in inFiles:
            binId = binIdFromFilename(fasta)
            
            binStats[binId] = {}
            seqStats[binId] = {}
            
            t = threading.Thread(target=self.__processBin, args=(fasta, outFolder, binStats[binId], seqStats[binId]))
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
        
        fout = open(os.path.join(outFolder, defaultValues.__CHECKM_DEFAULT_SEQ_STATS_FILE__), 'w')
        fout.write(str(seqStats))
        fout.close()
              
    def __processBin(self, fasta, outFolder, binStats, seqStats):
        """Thread safe processing of FASTA file."""
        logger = logging.getLogger()
        
        self.threadPool.acquire()
        try:
            self.varLock.acquire()
            try:
                self.numFilesStarted += 1
                self.logger.info('    Processing %s (%d of %d)' % (fasta, self.numFilesStarted, self.numFiles)) 
            finally:
                self.varLock.release()

            outDir = os.path.join(outFolder, os.path.basename(fasta))

            # read seqs
            seqs = readFasta(fasta)
            for seqId in seqs:
                seqStats[seqId] = {}
                
            # calculate GC statistics
            GC, stdGC = self.calculateGC(seqs, seqStats) 
            binStats['GC'] = GC
            binStats['GC std'] = stdGC
            
            # calculate statistics related to scaffold lengths
            maxScaffoldLen, maxContigLen, genomeSize, scaffold_N50, contig_N50, numContigs = self.calculateSeqLengthStats(seqs, seqStats)  
            binStats['# scaffolds'] = len(seqs)
            binStats['# contigs'] = numContigs
            binStats['Longest scaffold'] = maxScaffoldLen
            binStats['Longest contig'] = maxContigLen
            binStats['Genome size'] = genomeSize
            binStats['N50 (scaffolds)'] = scaffold_N50
            binStats['N50 (contigs)'] = contig_N50
            
            # calculate coding density statistics
            codingDensity, numORFs = self.calculateCodingDensity(outDir, genomeSize, seqStats)
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
            
    def calculateGC(self, seqs, seqStats):
        """Calculate fraction of nucleotides that are G or C."""
        totalGC = 0
        totalAT = 0
        gcPerSeq = []
        for seqId, seq in seqs.iteritems():  
            a, c, g, t = baseCount(seq)
            
            gc = g + c
            at = a + t
            
            totalGC += gc
            totalAT += at
            
            gcContent = float(gc) / (gc + at) 
            seqStats[seqId]['GC'] = gcContent
            
            if len(seq) > defaultValues.__CHECKM_DEFAULT_MIN_SEQ_LEN_GC_STD__:
                gcPerSeq.append(gcContent)
                
        GC = float(totalGC) / (totalGC + totalAT)
        
        varGC = mean(map(lambda x: (x - GC)**2, gcPerSeq))

        return GC, math.sqrt(varGC)
    
    def calculateSeqLengthStats(self, seqs, seqStats):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        scaffoldLens = []
        contigLens = []
        for scaffoldId, scaffold in seqs.iteritems():
            scaffoldLen = len(scaffold)
            scaffoldLens.append(scaffoldLen)
                        
            splitScaffold = scaffold.split(defaultValues.__CHECKM_DEFAULT_CONTIG_BREAK__)
            lenContigsInScaffold = []
            for contig in splitScaffold:
                contigLen = len(contig.replace('N', ''))
                if contigLen > 0:
                    lenContigsInScaffold.append(contigLen)

            contigLens += lenContigsInScaffold
            
            seqStats[scaffoldId]['Length'] = scaffoldLen
            seqStats[scaffoldId]['Total contig length'] = sum(lenContigsInScaffold)
            seqStats[scaffoldId]['# contigs'] = len(lenContigsInScaffold)
            
        scaffold_N50 = calculateN50(scaffoldLens)
        contig_N50 = calculateN50(contigLens)
            
        return max(scaffoldLens), max(contigLens), sum(scaffoldLens), scaffold_N50, contig_N50, len(contigLens)
        
    def calculateCodingDensity(self, outDir, genomeSize, seqStats):
        """Calculate coding density of putative genome bin."""
        ntFile = os.path.join(outDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_NT__)
        ntGenes = readFasta(ntFile)
        
        codingBasePairs = 0
        for geneId, gene in ntGenes.iteritems():
            codingBasePairs += len(gene)
            
            scaffoldId = geneId[0:geneId.rfind('_')]
            seqStats[scaffoldId]['# ORFs'] = seqStats[scaffoldId].get('# ORFs', 0) + 1
            seqStats[scaffoldId]['Coding bases'] = seqStats[scaffoldId].get('Coding bases', 0) + len(gene)
            
        return float(codingBasePairs) / genomeSize, len(ntGenes)
        
            