###############################################################################
#
# binTools.py - functions for exploring and modifying bins
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

import numpy as np

from common import binIdFromFilename, checkFileExists, readDistribution, findNearest
from seqUtils import readFasta, readFastaSeqIds, writeFasta, baseCount
from genomicSignatures import GenomicSignatures
from prodigal import ProdigalGeneFeatureParser

class BinTools():
    """Functions for exploring and modifying bins."""
    def __init__(self, threads=1):         
        self.logger = logging.getLogger()
    
    def __removeSeqs(self, seqs, seqsToRemove):
        """Remove sequences. """
        missingSeqIds = set(seqsToRemove).difference(set(seqs.keys()))
        if len(missingSeqIds) > 0:
            self.logger.error('  [Error] Missing sequence(s) specified for removal: ' + ', '.join(missingSeqIds) + '\n')
            sys.exit()
            
        for seqId in seqsToRemove:
            seqs.pop(seqId)
                
    def __addSeqs(self, seqs, refSeqs, seqsToAdd):
        """Add sequences. """
        missingSeqIds = set(seqsToAdd).difference(set(refSeqs.keys()))
        if len(missingSeqIds) > 0:
            self.logger.error('  [Error] Missing sequence(s) specified for addition: ' + ', '.join(missingSeqIds) + '\n')
            sys.exit()
        
        for seqId in seqsToAdd:
            seqs[seqId] = refSeqs[seqId]
    
    def modify(self, binFile, seqFile, seqsToAdd, seqsToRemove, outputFile):
        """Add and remove sequences from a file."""
        binSeqs = readFasta(binFile)
        
        # add sequences to bin
        if seqsToAdd != None:
            refSeqs = readFasta(seqFile)
            self.__addSeqs(binSeqs, refSeqs, seqsToAdd)
                
        # remove sequences from bin
        if seqsToRemove != None:
            self.__removeSeqs(binSeqs, seqsToRemove)
                
        # save modified bin
        writeFasta(binSeqs, outputFile)
        
    def removeOutliers(self, binFile, outlierFile, outputFile):
        """Remove sequences specified as outliers in the provided file."""
        
        binSeqs = readFasta(binFile)
        binIdToModify = binIdFromFilename(binFile)
        
        # get files to remove
        checkFileExists(outlierFile)
        seqsToRemove = []
        bHeader = True
        for line in open(outlierFile):
            if bHeader:
                bHeader = False
                continue
            
            lineSplit = line.split('\t')
            binId = lineSplit[0]
            
            if binId == binIdToModify:
                seqId = lineSplit[1]
                seqsToRemove.append(seqId)
                
        # remove sequences from bin
        if len(seqsToRemove) > 0:
            self.__removeSeqs(binSeqs, seqsToRemove)
        
        # save modified bin
        writeFasta(binSeqs, outputFile)

    def unique(self, binFiles):
        """Check if sequences are assigned to multiple bins."""
        
        # read seq ids from all bins
        binSeqs = {}
        for f in binFiles:
            binId = binIdFromFilename(f)
            binSeqs[binId] = readFastaSeqIds(f)

        # check for sequences assigned to multiple bins
        bDuplicates = False
        binIds = binSeqs.keys()
        for i in xrange(0, len(binIds)):
            for j in xrange(i+1, len(binIds)):
                seqInter = set(binSeqs[binIds[i]]).intersection(set(binSeqs[binIds[j]]))
                
                if len(seqInter) > 0:
                    bDuplicates = True
                    print '  Sequences shared between %s and %s: ' % (binIds[i], binIds[j]) + ', '.join(seqInter)
        
        if not bDuplicates:
            print '  No sequences assigned to multiple bins.'
            
    def gcDist(self, seqs):
        """GC statistics for bin."""
        GCs = []
        gcTotal = 0
        basesTotal = 0
        for _, seq in seqs.iteritems():
            a, c, g, t = baseCount(seq)
            gc = g + c
            bases = a + c + g + t
            
            GCs.append(float(gc) / (bases))
            
            gcTotal += gc
            basesTotal += bases

        meanGC = float(gcTotal) / basesTotal
        deltaGCs = np.array(GCs) - meanGC
        
        return meanGC, deltaGCs, GCs
    
    def codingDensityDist(self, seqs, prodigalParser):
        """Coding density statistics for bin."""
        CDs = []
        seqLens = []
        
        codingBasesTotal = 0
        basesTotal = 0
        for seqId, seq in seqs.iteritems():
            seqLens.append(len(seq))

            codingBases = prodigalParser.codingBases(seqId)
                
            a, c, g, t = baseCount(seq)
            validBases = a + c +g + t
            CDs.append(float(codingBases) / validBases)
            
            codingBasesTotal += codingBases
            basesTotal += validBases
            
            seqLens.append(len(seq))

        meanCD = float(codingBasesTotal) / basesTotal
        deltaCDs = np.array(CDs) - meanCD
        
        return meanCD, deltaCDs, CDs
    
    def binTetraSig(self, seqs, tetraSigs):
        """Tetranucleotide signature for bin. """
        binSize = 0
        for _, seq in seqs.iteritems():
            binSize += len(seq)
        
        binSig = None
        for seqId, seq in seqs.iteritems():
            weightedTetraSig = tetraSigs[seqId] * (float(len(seq)) / binSize)
            if binSig == None:
                binSig = weightedTetraSig
            else:
                binSig += weightedTetraSig
                
        return binSig
    
    def tetraDiffDist(self, seqs, genomicSig, tetraSigs, binSig):
        """TD statistics for bin."""
        deltaTDs = np.zeros(len(seqs))
        for i, seqId in enumerate(seqs.keys()):
            deltaTDs[i] = genomicSig.distance(tetraSigs[seqId], binSig)
        
        return np.mean(deltaTDs), deltaTDs
            
    def identifyOutliers(self, outFolder, binFiles, tetraProfileFile, distribution, reportType, outputFile):
        """Identify sequences that are outliers."""   
        gcBounds = readDistribution(distribution, 'gc_dist')
        cdBounds = readDistribution(distribution, 'cd_dist')
        tdBounds = readDistribution(distribution, 'td_dist')

        fout = open(outputFile, 'w')
        fout.write('Bin Id\tSequence Id\tSequence length\tOutlying distributions')
        fout.write('\tSequence GC\tMean bin GC\tLower GC bound (%s%%)\tUpper GC bound (%s%%)' % (distribution, distribution))
        fout.write('\tSequence CD\tMean bin CD\tLower CD bound (%s%%)\tUpper CD bound (%s%%)' % (distribution, distribution))
        fout.write('\tSequence TD\tMean bin TD\tUpper TD bound (%s%%)\n' % distribution)
        
        processedBins = 0
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)
            
            processedBins += 1
            self.logger.info('  Finding outliers in bin %s (%d of %d)' % (binId, processedBins, len(binFiles)))
                
            seqs = readFasta(binFile)
            
            meanGC, deltaGCs, seqGC = self.gcDist(seqs)
            
            genomicSig = GenomicSignatures(K=4, threads=1)
            tetraSigs = genomicSig.read(tetraProfileFile)
            binSig = self.binTetraSig(seqs, tetraSigs)
            meanTD, deltaTDs = self.tetraDiffDist(seqs, genomicSig, tetraSigs, binSig)
            
            gffFile = os.path.join(outFolder, os.path.basename(binFile), 'prodigal.gff')
            prodigalParser = ProdigalGeneFeatureParser(gffFile)
            meanCD, deltaCDs, CDs = self.codingDensityDist(seqs, prodigalParser)
            
            index = 0
            for seqId, seq in seqs.iteritems(): 
                closestGC = findNearest(np.array(gcBounds.keys()), meanGC)
                closestWindowSize = findNearest(np.array(gcBounds[closestGC].keys()), len(seq))
                gcLowerBound = gcBounds[closestGC][closestWindowSize][0]
                gcUpperBound = gcBounds[closestGC][closestWindowSize][1]
                
                closestCD = findNearest(np.array(cdBounds.keys()), meanCD)
                closestWindowSize = findNearest(np.array(cdBounds[closestCD].keys()), len(seq))
                cdLowerBound = cdBounds[closestCD][closestWindowSize][0]
                cdUpperBound = cdBounds[closestCD][closestWindowSize][1]
     
                closestWindowSize = findNearest(np.array(tdBounds.keys()), len(seq))
                tdBound = tdBounds[closestWindowSize]

                outlyingDists = []
                if deltaGCs[index] < gcLowerBound or deltaGCs[index] > gcUpperBound:
                    outlyingDists.append('GC')
                    
                if deltaCDs[index] < cdLowerBound or deltaCDs[index] > cdUpperBound:
                    outlyingDists.append('CD')
                    
                if deltaTDs[index] > tdBound:
                    outlyingDists.append('TD')
                     
                if (reportType == 'any' and len(outlyingDists) >= 1) or (reportType == 'all' and len(outlyingDists) == 3):
                    fout.write(binId + '\t' + seqId + '\t%d' % len(seq) + '\t' + ','.join(outlyingDists))
                    fout.write('\t%.1f\t%.1f\t%.1f\t%.1f' % (seqGC[index]*100, meanGC*100, (meanGC+gcLowerBound)*100, (meanGC+gcUpperBound)*100))
                    fout.write('\t%.1f\t%.1f\t%.1f\t%.1f' % (CDs[index]*100, meanCD*100, (meanCD+cdLowerBound)*100, (min(1.0, meanCD+cdUpperBound))*100))
                    fout.write('\t%.1f\t%.1f\t%.1f' % (deltaTDs[index]*100, meanTD*100, tdBound*100) + '\n')
                    
                index += 1
                              
        fout.close()
                    