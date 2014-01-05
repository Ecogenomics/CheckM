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
import json
import multiprocessing as mp

import numpy as np

class CalculateBounds(object):
    def __init__(self):
        self.lock = mp.Lock()
        self.counter = mp.Manager().Value('i', 0)
        
    def run(self, genomeDir, outputFile, stepSize, width, minGenomes, threads):
        """Process central values in parallel."""  
        
        # get mean value of all genomes
        print 'Parsing mean values from all genomes.'
        files = os.listdir(genomeDir)
        
        meanValues = {}
        for f in files:
            if not f.endswith('.tsv'):
                continue
            
            with open(os.path.join(genomeDir, f), 'r') as fin:
                dataLine = fin.readline()
                
                genomeId = f[0:f.rfind('.')]
                meanValue = float(dataLine.split('=')[1])
                meanValues[genomeId] = meanValue
                  
        # process different mean values in parallel   
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        centreValues = np.arange(0.0, 1.0 + 0.5*stepSize, stepSize)
        for v in centreValues:
            workerQueue.put(v)

        for _ in range(threads):
            workerQueue.put(None)

        CIs = np.arange(0, 100 + 0.5, 0.5).tolist()
        calcProc = [mp.Process(target = self.__calculateResults, args = (genomeDir, CIs, meanValues, width, minGenomes, workerQueue, writerQueue)) for _ in range(threads)]
        writeProc = mp.Process(target = self.__storeResults, args = (outputFile, len(centreValues), writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None))
        writeProc.join()
        
    def __getGenomesInRange(self, meanValues, minValue, maxValue):
        genomeIds = []
        for genomeId, value in meanValues.iteritems():
            if value >= minValue and value <= maxValue:
                genomeIds.append(genomeId)
                
        return genomeIds
    
    def __getDeltaForWindows(self, genomeDir, genomeIds):
        d = {}
        for genomeId in genomeIds:
            bReadDist = False
            for line in open(os.path.join(genomeDir, genomeId + '.tsv')):
                if 'Windows Size' in line:
                    windowSize = int(line.split('=')[1].strip())
                    if windowSize not in d:
                        d[windowSize] = []
                    bReadDist = True
                elif bReadDist:
                    bReadDist = False
                    d[windowSize].extend([float(x) for x in line.split(',')])
                
        return d
            
    def __calculateResults(self, genomeDir, CIs, meanValues, width, minGenomes, queueIn, queueOut):
        """Calculate confidence interval for a specific mean value (e.g., GC) and different sequence lengths."""
        while True:
            meanValue = queueIn.get(block=True, timeout=None) 
            if meanValue == None:
                break

            # get genomes within a given width of mean value
            genomeIds = self.__getGenomesInRange(meanValues, meanValue-width, meanValue+width)
            
            if len(genomeIds) < minGenomes:
                with self.lock:
                    self.counter.value += 1
                continue
            
            windows = self.__getDeltaForWindows(genomeDir, genomeIds)  
            d = {} # Note: defaultdict doesn't play well with 'ast' library
            for windowSize, testPts in windows.iteritems():
                testPts = np.array(testPts)
                
                d[windowSize] = {}
                percentiles = np.percentile(testPts, CIs)
                for ci, p in zip(CIs, percentiles):
                    d[windowSize][ci] = p
                
            queueOut.put((meanValue, d))
                    
            with self.lock:
                self.counter.value += 1

    def __storeResults(self, outputFile, totalItems, queue):
        """Store confidence intervals (i.e., to shared memory)."""
        print 'Calculating distributions for %d mean values.' % totalItems
        
        globalDist = {} # Note: defaultdict doesn't play well with 'ast' library
        while True:
            meanValue, dist = queue.get(block=True, timeout=None)
            if meanValue == None:
                break

            statusStr = '  Finished processing %d of %d (%.2f%%) steps.' % (self.counter.value, totalItems, float(self.counter.value)*100/totalItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            globalDist[meanValue] = dist 

        sys.stdout.write('\n')
        
        fout = open(outputFile, 'w')
        fout.write(str(globalDist))
        fout.close()
              
if __name__ == '__main__':  
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                            description='Calculate distribution bounds.')
    parser.add_argument('genome_dir', help='directory with distributions for each genome')
    parser.add_argument('output_file', help='output file for distributions')
    parser.add_argument('-s', '--step_size', help='step size for calculating distribution', type = float, default = 0.01)
    parser.add_argument('-w', '--width', help='width around mean of distribution', type = float, default = 0.015)
    parser.add_argument('--min_genomes', help='minimum genomes to calculate statistics', type = int, default = 5)
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)
    
    args = parser.parse_args()
    
    calculateBounds = CalculateBounds()
    calculateBounds.run(args.genome_dir, args.output_file, args.step_size, args.width, args.min_genomes, args.threads)