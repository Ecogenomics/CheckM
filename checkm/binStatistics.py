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
import multiprocessing as mp
import math
import logging

import defaultValues
from seqUtils import readFasta, baseCount, calculateN50
from common import binIdFromFilename, makeSurePathExists

from numpy import mean

class BinStatistics():
    """Calculate statistics (GC, coding density, etc.) for genome bins."""
    def __init__(self, threads):   
        self.logger = logging.getLogger()        
        self.totalThreads = threads

    def calculate(self, binFiles, outDir, binStatsFile, seqStatsFile):
        """Calculate statistics for each putative genome bin."""
        
        # process each bin
        self.logger.info("  Calculating genome statistics for %d bins with %d threads:" % (len(binFiles), self.totalThreads))

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for binFile in binFiles:
            workerQueue.put(binFile)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__processBin, args = (outDir, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportProgress, args = (outDir, binStatsFile, seqStatsFile, len(binFiles), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put((None, None, None))
        writeProc.join()
              
    def __processBin(self, outDir, queueIn, queueOut):
        """Thread safe bin processing."""
        while True:    
            binFile = queueIn.get(block=True, timeout=None) 
            if binFile == None:
                break   
            
            binStats = {}
            scaffoldStats = {}
            
            binId = binIdFromFilename(binFile)
            binDir = os.path.join(outDir, 'bins', binId)
            makeSurePathExists(binDir)

            # read scaffolds
            scaffolds = readFasta(binFile) 
            for seqId in scaffolds:
                scaffoldStats[seqId] = {}
                
            # calculate GC statistics
            GC, stdGC = self.calculateGC(scaffolds, scaffoldStats) 
            binStats['GC'] = GC
            binStats['GC std'] = stdGC
            
            # calculate statistics related to scaffold lengths
            maxScaffoldLen, maxContigLen, genomeSize, scaffold_N50, contig_N50, numContigs = self.calculateScaffoldLengthStats(scaffolds, scaffoldStats)  
            binStats['# scaffolds'] = len(scaffolds)
            binStats['# contigs'] = numContigs
            binStats['Longest scaffold'] = maxScaffoldLen
            binStats['Longest contig'] = maxContigLen
            binStats['Genome size'] = genomeSize
            binStats['N50 (scaffolds)'] = scaffold_N50
            binStats['N50 (contigs)'] = contig_N50
            
            # calculate coding density statistics            
            codingDensity, numORFs = self.calculateCodingDensity(binDir, genomeSize, scaffoldStats)
            binStats['Coding density'] = codingDensity
            binStats['# predicted ORFs'] = numORFs
            
            queueOut.put((binId, binStats, scaffoldStats))
            
    def __reportProgress(self, outDir, binStatsFile, seqStatsFile, numBins, queueIn):
        """Report number of processed bins and write statistics to file."""      
        
        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        
        binStats = {}
        seqStats = {}
        while True:
            binId, curBinStats, curSeqStats = queueIn.get(block=True, timeout=None)
            if binId == None:
                break
            
            binStats[binId] = curBinStats
            seqStats[binId] = curSeqStats
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
                
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
       
        # save results
        storagePath = os.path.join(outDir, 'storage')  
        
        fout = open(os.path.join(storagePath, binStatsFile), 'w')
        fout.write(str(binStats))
        fout.close()
        
        fout = open(os.path.join(storagePath, seqStatsFile), 'w')
        fout.write(str(seqStats))
        fout.close()
         
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
            
            if len(seq) > defaultValues.MIN_SEQ_LEN_GC_STD:
                gcPerSeq.append(gcContent)
                
        GC = float(totalGC) / (totalGC + totalAT)
        
        varGC = 0
        if len(gcPerSeq) > 1:
            varGC = mean(map(lambda x: (x - GC)**2, gcPerSeq))

        return GC, math.sqrt(varGC)
    
    def calculateScaffoldLengthStats(self, scaffolds, scaffoldsStats):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        scaffoldLens = []
        contigLens = []
        for scaffoldId, scaffold in scaffolds.iteritems():
            scaffoldLen = len(scaffold)
            scaffoldLens.append(scaffoldLen)
                        
            splitScaffold = scaffold.split(defaultValues.CONTIG_BREAK)
            lenContigsInScaffold = []
            for contig in splitScaffold:
                contigLen = len(contig.replace('N', ''))
                if contigLen > 0:
                    lenContigsInScaffold.append(contigLen)

            contigLens += lenContigsInScaffold
            
            scaffoldsStats[scaffoldId]['Length'] = scaffoldLen
            scaffoldsStats[scaffoldId]['Total contig length'] = sum(lenContigsInScaffold)
            scaffoldsStats[scaffoldId]['# contigs'] = len(lenContigsInScaffold)
            
        scaffold_N50 = calculateN50(scaffoldLens)
        contig_N50 = calculateN50(contigLens)
            
        return max(scaffoldLens), max(contigLens), sum(scaffoldLens), scaffold_N50, contig_N50, len(contigLens)
        
    def calculateCodingDensity(self, outDir, genomeSize, seqStats):
        """Calculate coding density of putative genome bin."""
        ntFile = os.path.join(outDir, defaultValues.PRODIGAL_NT)
        ntGenes = readFasta(ntFile)
        
        codingBasePairs = 0
        for geneId, gene in ntGenes.iteritems():
            codingBasePairs += len(gene)
            
            scaffoldId = geneId[0:geneId.rfind('_')]
            seqStats[scaffoldId]['# ORFs'] = seqStats[scaffoldId].get('# ORFs', 0) + 1
            seqStats[scaffoldId]['Coding bases'] = seqStats[scaffoldId].get('Coding bases', 0) + len(gene)
            
        return float(codingBasePairs) / genomeSize, len(ntGenes)
        
            