#!/usr/bin/env python

###############################################################################
#
# calculateBoundsTD.py - find confidence intervals for tetranucleotide
#                            signature differences
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
import argparse
import multiprocessing as mp

class CalculateBoundsTD(object):
    def __init__(self):
        self.manager = mp.Manager()
        self.dist = self.manager.dict()
        
        self.lock = mp.Lock()
        self.counter = self.manager.Value('i', 0)
        
    def run(self, genomeDir, outputDir):
        
        # get all genome ids
        print 'Determine genome IDs.'
        files = os.listdir(genomeDir)
        genomeIds = []
        for f in files:
            genomeIds.append(f[0:f.rfind('.')])
            
        print '  Total genomes: ' + str(len(genomeIds))

        # get tetranucleotide signature differences for each window size
        print 'Reading windows.'
        windows = {}
        processedGenomes = 0
        for genomeId in genomeIds:
            processedGenomes += 1
            statusStr = '  Finished processing %d of %d (%.2f%%) genomes' % (processedGenomes, len(genomeIds), float(processedGenomes)*100/len(genomeIds))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            bReadDist = False
            for line in open(os.path.join(genomeDir, genomeId + '.txt')):
                if 'Windows Size' in line:
                    windowSize = int(line.split('=')[1].strip())
                    bReadDist = True
                elif bReadDist:
                    bReadDist = False
                    windows[windowSize] = windows.get(windowSize, []) + [float(x) for x in line.split(',')]
        sys.stdout.write('\n')
                
        # get CI for each window size
        print 'Determining CI for each window size.'
        CIs = [0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0]
        dist = {}
        processedWindows = 0
        for windowSize in windows:    
            processedWindows += 1
            statusStr = '  Finished processing %d of %d (%.2f%%) window sizes; number of windows = %d' % (processedWindows, len(windows), float(processedWindows)*100/len(windows), len(windows[windowSize]))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            sortedWindows = sorted(windows[windowSize])
            
            for ci in CIs:
                if ci == 0:
                    upperIndex = len(sortedWindows) - 1
                else:
                    upperIndex = int((1.0 - ci) * len(sortedWindows) + 0.5) - 1
                          
                windowDict = dist.get(ci, {})
                windowDict[windowSize] = sortedWindows[upperIndex]
                dist[ci] = windowDict          
        sys.stdout.write('\n')               
                    
        # save results to file
        print 'Writing results to file.'
        for ci in CIs:
            ciStr = '%.2f' % ((1.0 - ci)*100)
            ciStr = ciStr.rstrip('0').rstrip('.')
            fout = open(outputDir + '/distribution_%s.txt' % ciStr, 'w')
            fout.write(str(dist[ci]))
            fout.close()
            
if __name__ == '__main__':  
    parser = argparse.ArgumentParser(description='Calculate distribution bounds.')
    parser.add_argument('genome_dir', help='directory with distributions for each genome')
    parser.add_argument('output_dir', help='directory to write out distributions')
    
    args = parser.parse_args()
    
    calculateBoundsTD = CalculateBoundsTD()
    calculateBoundsTD.run(args.genome_dir, args.output_dir)