#!/usr/bin/env python

###############################################################################
#
# calculateBounds.py - find confidence intervals for distributions
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

import numpy as np

class CalculateBounds(object):
    def __init__(self):
        self.manager = mp.Manager()
        self.dist = self.manager.dict()
        
        self.lock = mp.Lock()
        self.counter = self.manager.Value('i', 0)
        
    def run(self, genomeDir, outputDir, stepSize, width, minGenomes, threads):
        """Process central values in parallel."""  
        
        # get mean value of all genomes
        print 'Parsing mean values from all genomes.'
        files = os.listdir(genomeDir)
        
        meanValues = {}
        for f in files:
            if not f.endswith('.tsv'):
                continue
            
            with open(os.path.join(genomeDir, f), 'r') as fin:
                gcLine = fin.readline()
                
                genomeId = f[0:f.rfind('.')]
                meanValue = float(gcLine.split('=')[1])
                meanValues[genomeId] = meanValue
                  
        # process different mean values in parallel   
        CIs = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.0005, 0]
        centreValues = np.arange(0.0, 1.0 + 0.5*stepSize, stepSize)

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for v in centreValues:
            workerQueue.put(v)

        for _ in range(threads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__calculateResults, args = (genomeDir, CIs, meanValues, width, minGenomes, workerQueue, writerQueue)) for _ in range(threads)]
        writeProc = mp.Process(target = self.__storeResults, args = (writerQueue, len(centreValues)))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None, None, None, None))
        writeProc.join()
        
        for ci in CIs:
            ciStr = '%.2f' % ((1.0 - 2*ci)*100)
            ciStr = ciStr.rstrip('0').rstrip('.')
            fout = open(outputDir + '/distribution_%s.txt' % ciStr, 'w')
            fout.write(str(self.dist[ci]))
            fout.close()
        
    def __getGenomesInRange(self, meanValues, minValue, maxValue):
        genomeIds = []
        for genomeId, gc in meanValues.iteritems():
            if gc >= minValue and gc <= maxValue:
                genomeIds.append(genomeId)
                
        return genomeIds
    
    def __getDeltaForWindows(self, genomeDir, genomeIds):
        d = {}
        for genomeId in genomeIds:
            bReadDist = False
            for line in open(os.path.join(genomeDir, genomeId + '.tsv')):
                if 'Windows Size' in line:
                    windowSize = int(line.split('=')[1].strip())
                    bReadDist = True
                elif bReadDist:
                    bReadDist = False
                    d[windowSize] = d.get(windowSize, []) + [float(x) for x in line.split(',')]
                
        return d
            
    def __calculateResults(self, genomeDir, CIs, meanValues, width, minGenomes, queueIn, queueOut):
        """Calculate confidence interval for a specific mean GC and different sequence lengths."""
        while True:
            meanValue = queueIn.get(block=True, timeout=None) 
            if meanValue == None:
                break

            # get genomes within GC width of mean GC
            genomeIds = self.__getGenomesInRange(meanValues, meanValue-width, meanValue+width)
            
            if len(genomeIds) < minGenomes:
                with self.lock:
                    self.counter.value += 1
                continue
            
            windows = self.__getDeltaForWindows(genomeDir, genomeIds)     
            for windowSize in windows:
                sortedWindows = sorted(windows[windowSize])
                
                for ci in CIs:
                    if ci == 0:
                        lowerIndex = 0
                        upperIndex = len(sortedWindows) - 1
                    else:
                        lowerIndex = int(ci * len(sortedWindows)) - 1
                        upperIndex = int((1.0 - ci) * len(sortedWindows) + 0.5) - 1
                
                    queueOut.put((ci, meanValue, windowSize, sortedWindows[lowerIndex], sortedWindows[upperIndex]))
                    
            with self.lock:
                self.counter.value += 1

    def __storeResults(self, queue, totalItems):
        """Store confidence intervals (i.e., to shared memory)."""
        print 'Calculating distributions for %d mean values.' % totalItems
        while True:
            ci, meanValue, windowSize, lowerGC, upperGC = queue.get(block=True, timeout=None)
            if ci == None:
                break

            statusStr = '  Finished processing %d of %d (%.2f%%) steps; last value = %.2f' % (self.counter.value, totalItems, float(self.counter.value)*100/totalItems, meanValue)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            # this pattern is necessary for the Manager to serialize data
            valueDict = self.dist.get(ci, {})
            windowDict = valueDict.get(meanValue, {})
            windowDict[windowSize] = [lowerGC, upperGC]
            valueDict[meanValue] = windowDict
            self.dist[ci] = valueDict       

        sys.stdout.write('\n')
              
if __name__ == '__main__':  
    parser = argparse.ArgumentParser(description='Calculate distribution bounds.')
    parser.add_argument('genome_dir', help='directory with distributions for each genome')
    parser.add_argument('output_dir', help='directory to write out distributions')
    parser.add_argument('-s', '--step_size', help='step size for calculating distribution', type = float, default = 0.01)
    parser.add_argument('-w', '--width', help='width around mean of distribution', type = float, default = 0.015)
    parser.add_argument('--min_genomes', help='minimum genomes to calculate delta GC statistics', type = int, default = 5)
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)
    
    args = parser.parse_args()
    
    calculateBounds = CalculateBounds()
    calculateBounds.run(args.genome_dir, args.output_dir, args.step_size, args.width, args.min_genomes, args.threads)