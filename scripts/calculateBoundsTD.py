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
import random
import argparse
import pickle

import numpy as np

class CalculateBoundsTD(object):
    def __init__(self):
        pass
        
    def run(self, genomeDir, outputFile):
        
        # get all genome ids
        print 'Determine genome IDs.'
        files = os.listdir(genomeDir)
        genomeIds = []
        
        # for unknown reasons these genomes have NaN values so
        # are being ignored
        for f in files:
            genomeId = f[0:f.rfind('.')]
            genomeIds.append(genomeId)
            
        print '  Total genomes: ' + str(len(genomeIds))

        # get tetranucleotide signature differences for each window size
        print 'Reading windows.'
        windows = {} # Note: defaultdict doesn't play well with 'ast' library
        processedGenomes = 0
        badGenomes = set()
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
                    
                    # minimize influence of any given genome and 
                    # keep number of samples to a reasonable level
                    testPts = [float(x) for x in line.split(',')]
                    
                    # there are a few genomes with NaN's in there data file
                    # (these are currently just being ignored)
                    if np.isnan(np.sum(testPts)):
                        badGenomes.add(genomeId)
                        continue
                        
                    rndTestPts = random.sample(testPts, min(10000, testPts)) 
                    
                    if windowSize not in windows:
                        windows[windowSize] = []
                    windows[windowSize].extend(rndTestPts)
        sys.stdout.write('\n')
        
        print ''
        print '# bad genomes: ' + str(len(badGenomes))
        print badGenomes
        print ''
                
        # get CI for each window size
        print 'Determining CI for each window size.'
        CIs = np.arange(0, 100 + 0.5, 0.5).tolist()
        d = {} # Note: defaultdict doesn't play well with 'ast' library
        processedWindows = 0
        for windowSize, testPts in windows.iteritems():  
            processedWindows += 1
            statusStr = '  Finished processing %d of %d (%.2f%%) window sizes; number of windows = %d' % (processedWindows, len(windows), float(processedWindows)*100/len(windows), len(windows[windowSize]))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            testPts = np.array(testPts)
            
            # where are the nan's coming from??
            if np.isnan(np.sum(rndTestPts)):
                print 'Sneaky NaN: ' + windowSize
                sys.exit()
            
            d[windowSize] = {}
            percentiles = np.percentile(testPts, CIs)
            
            # where are the nan's coming from??
            if np.isnan(np.sum(percentiles)):
                print 'NaN in percentiles: ' + windowSize
                print len(testPts)
                print percentiles
                sys.exit()
            
            for ci, p in zip(CIs, percentiles):
                d[windowSize][ci] = p
                                
        sys.stdout.write('\n')               
                    
        # save results to file
        print 'Writing results to file.'
        fout = open(outputFile, 'w')
        fout.write(str(d))
        fout.close()
            
if __name__ == '__main__':  
    parser = argparse.ArgumentParser(description='Calculate distribution bounds.')
    parser.add_argument('genome_dir', help='directory with distributions for each genome')
    parser.add_argument('output_file', help='output file for distributions')
    
    args = parser.parse_args()
    
    calculateBoundsTD = CalculateBoundsTD()
    calculateBoundsTD.run(args.genome_dir, args.output_file)