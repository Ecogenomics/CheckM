#!/usr/bin/env python

###############################################################################
#
# distributionDeltaGC.py - delta GC distribution over reference genomes
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
import argparse
import string
from random import randint

import numpy as np

from checkm.seqUtils import readFasta, baseCount

class DeltaGC(object):
    def __init__(self):
        pass
    
    def getGenomesInGC_Range(self, genomeGC, minGC, maxGC):
        genomeIds = []
        for genomeId, gc in genomeGC.iteritems():
            if gc >= minGC and gc < maxGC:
                genomeIds.append(genomeId)
                
        return genomeIds
    
    def calculateGC(self, genomeId, genomeDir):
        # calculate GC of genome and make sure it corresponds to metadata
        seqs = readFasta(genomeDir + '/' + genomeId + '/' + genomeId + '.fna')
        
        gc = 0
        totalBases = 0
        for seq in seqs.values():
            a, c, g, t = baseCount(seq)
            gc += (g + c)
            totalBases += (a+c+g+t)
        
        return float(gc) / totalBases
    
    def run(self, metadataFile, genomeDir, gcStepSize, gcWidth, minLen, maxLen, lenStepSize, numWindows, minGenomes):
        # read metadata file
        print 'Determining finished prokaryotic reference genomes.'
        genomeGC = {}
        bHeader = True
        for line in open(metadataFile):
            lineSplit = [x.strip() for x in line.split('\t')]
            
            if bHeader:
                bHeader = False
                continue
            
            domain = lineSplit[1]
            status = lineSplit[2]
            if status == 'Finished' and (domain == 'Bacteria' or domain == 'Archaea'):
                genomeId = lineSplit[0]
                gc = float(lineSplit[69])
                
                if os.path.exists(genomeDir + '/' + genomeId + '/' + genomeId + '.fna'):
                    genomeGC[genomeId] = gc
                    
        print '  Identified reference genomes: ' + str(len(genomeGC))

        # sanity check calculate GC for each genome
        if False:
            print 'Calculating GC of each genome.'
            for genomeId in genomeGC:
                gc = self.calculateGC(genomeId, genomeDir)
                
                print gc, genomeGC[genomeId]
                if gc < genomeGC[genomeId] - 0.01 or gc > genomeGC[genomeId] + 0.01:
                    sys.stderr.write('[Error] GC in metadata is incorrect from genome: ' + genomeId)
                    sys.exit()
        
        # calculate delta-GC distribution for genomes with different mean GC and
        # over varying window sizes
        print '\nCalculating delta-GC statistics:'
        for gc in np.arange(0.0, 1.0 + 0.5*gcStepSize, gcStepSize):
            print '  GC = ' + str(gc)
            
            # get genomes for each GC step
            genomeIds = self.getGenomesInGC_Range(genomeGC, gc-gcWidth, gc+gcWidth)
            print '  Genomes in range %.2f +/- %.3f: %d' % (gc, gcWidth, len(genomeIds))
            
            if len(genomeIds) > minGenomes:
                print '  Passed minimum genome check.'
            else:
                print '  Failed minimum genome check.'
                print '  ----------------------------------------------------------'
                continue
            
            fout = open('./deltaGC/gc_' + str(gc) + '.txt', 'w')
            fout.write('Mean GC = ' + str(gc) + '\n')
            fout.write('Genomes = ' + str(len(genomeIds)) + '\n')
            fout.write('GC Width = ' + str(gcWidth) + '\n')
            fout.write('# windows = ' + str(numWindows) + '\n\n')
            
            print '  Calculating delta-GC statistics.'
            numProcessedGenomes = 0
            for genomeId in genomeIds:
                fout.write('Genome ID = ' + genomeId + '\n')
                numProcessedGenomes += 1
                processedGenomesStr = '    Processing genome %d of %d' % (numProcessedGenomes, len(genomeIds))
                sys.stdout.write('%s\r' % processedGenomesStr)
                sys.stdout.flush()
                
                seqs = readFasta(genomeDir + '/' + genomeId + '/' + genomeId + '.fna')
                
                # make single seqs
                genomeSeq = ''.join(seqs.values()).upper()
                
                # create numpy array of genome with G and C as 0, A and T as 1, and degenerate/ambiguous bases as 2
                trans = string.maketrans('CGATUNRYMKBDHVWS', '0011122222222222')
                transGenomeSeq = genomeSeq.translate(trans)
                genomeList = np.array(map(int, transGenomeSeq), dtype=np.int)
                  
                meanGC = genomeGC[genomeId]
                for windowSize in xrange(minLen, maxLen + 1, lenStepSize):
                    endWindowPos = len(genomeSeq) - windowSize
                    requiredBasePairs = 0.9*windowSize
                    
                    deltaGCs = []
                    while len(deltaGCs) != numWindows:
                        # pick random window
                        startWindow = randint(0, endWindowPos)
                        
                        # calculate GC
                        counts = np.bincount(genomeList[startWindow:(startWindow+windowSize)])
                        gc = counts[0]
                        totalBases = gc + counts[1]
                        
                        if totalBases < requiredBasePairs:
                            # there are N's in the window so skip it
                            continue
                        
                        gcPer = float(gc) / totalBases
                        deltaGCs.append(meanGC - gcPer)
                        
                    fout.write('Windows Size = ' + str(windowSize) + '\n')
                    fout.write(','.join(map(str, deltaGCs)) + '\n')
            fout.close()
            print '\n  ----------------------------------------------------------'
                    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate delta GC distribution over reference genomes.')
    parser.add_argument('metadata_file', help='IMG metadata file.')
    parser.add_argument('genome_dir', help='IMG genome directory.')
    parser.add_argument('-g', '--gc_step_size', help='GC step size', type = float, default = 0.01)
    parser.add_argument('-w', '--gc_width', help='GC bin width', type = float, default = 0.02)
    parser.add_argument('--min_len', help='Minimum fragment length', type = int, default = 500)
    parser.add_argument('--max_len', help='Maximum fragment length', type = int, default = 10000)
    parser.add_argument('--len_step_size', help='Length step size', type = int, default = 500)
    parser.add_argument('--num_windows', help='Number of windows to sample', type = int, default = 10000)
    parser.add_argument('--min_genomes', help='Minimum genomes to calculate delta GC statistics', type = int, default = 5)
    
    args = parser.parse_args()
    
    deltaGC = DeltaGC()
    deltaGC.run(args.metadata_file, args.genome_dir, args.gc_step_size, args.gc_width, args.min_len, args.max_len, args.len_step_size, args.num_windows, args.min_genomes)