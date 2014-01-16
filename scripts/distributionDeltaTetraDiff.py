#!/usr/bin/env python

###############################################################################
#
# distributionDeltaTetra.py - calculate distribution of tetranucleotide
#                             distances for all reference genomes
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

import time
import sys
import os
import argparse
from random import randint
import multiprocessing as mp

import numpy as np

from checkm.seqUtils import readGenomicSeqsFromFasta
from checkm.genomicSignatures import GenomicSignatures

class DeltaTetraDiff(object):
    def __init__(self):
        pass

    def __calculateResults(self, windowSizes, numWindows, genomeDir, queueIn, queueOut):
        while True:
            genomeId = queueIn.get(block=True, timeout=None)
            if genomeId == None:
                break

            start = time.time()

            seqs = readGenomicSeqsFromFasta(os.path.join(genomeDir, genomeId, genomeId + '.fna'))
            genomeScaffold = 'NNNN'.join(seqs.values())

            # calculate tetranucleotide signature of genome
            gsCalculator = GenomicSignatures(4)
            genomeSig = gsCalculator.seqSignature(genomeScaffold)

            fout = open('./deltaTD/' + genomeId + '.tsv', 'w')
            fout.write('# Tetra signature = ' + str(genomeSig) + '\n')
            fout.close()
            sys.exit()

            # calculate tetranucleotide distance distribution for different window sizes
            startW = time.time()
            for windowSize in windowSizes:
                endWindowPos = len(genomeScaffold) - windowSize
                if endWindowPos <= 0:
                    # This might occur for the largest window sizes and smallest genomes
                    break

                deltaTDs = []
                while len(deltaTDs) != numWindows:
                    # pick random window
                    startWindow = randint(0, endWindowPos)

                    windowSig = gsCalculator.seqSignature(genomeScaffold[startWindow:(startWindow+windowSize)])
                    dist = gsCalculator.distance(genomeSig, windowSig)
                    deltaTDs.append(dist)

                fout.write('Windows Size = ' + str(windowSize) + '\n')
                fout.write(','.join(map(str, deltaTDs)) + '\n')
            fout.close()
            endW = time.time()
            print endW - startW

            queueOut.put(genomeId)

            end = time.time()
            print end - start

    def __storeResults(self, queue, numGenomes):
        processedRef = 0
        while True:
            genomeId = queue.get(block=True, timeout=None)
            if genomeId == None:
                break

            processedRef += 1
            processedStr = 'Finished processing %d of %d (%.2f%%) genomes.' % (processedRef, numGenomes, float(processedRef)*100/numGenomes)
            sys.stdout.write('%s\r' % processedStr)
            sys.stdout.flush()
        sys.stdout.write('\n')

    def run(self, metadataFile, genomeDir, numWindows, numThreads):
        # read metadata file
        print 'Determining finished prokaryotic reference genomes.'
        genomeIds = []
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

                if os.path.exists(os.path.join(genomeDir, genomeId, genomeId + '.fna')):
                    genomeIds.append(genomeId)

        genomeIds = genomeIds[0:1]

        print '  Identified reference genomes: ' + str(len(genomeIds))

        # window sizes to sample
        windowSizes = [ws for ws in np.arange(500, 1000, 100)]
        windowSizes += [ws for ws in np.arange(1000, 2000, 200)]
        windowSizes += [ws for ws in np.arange(2000, 5000, 500)]
        windowSizes += [ws for ws in np.arange(5000, 10000, 1000)]
        windowSizes += [ws for ws in np.arange(10000, 50000, 5000)]
        windowSizes += [ws for ws in np.arange(50000, 100000, 10000)]
        windowSizes += [ws for ws in np.arange(100000, 400000, 100000)]
        windowSizes += [ws for ws in np.arange(400000, 1000001, 200000)]

        print '# window sizes: ' + str(len(windowSizes))

        # sample windows from each genome
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for genomeId in genomeIds:
            workerQueue.put(genomeId)

        for _ in range(numThreads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__calculateResults, args = (windowSizes, numWindows, genomeDir, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target = self.__storeResults, args = (writerQueue, len(genomeIds)))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put(None)
        writeProc.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate GC distribution over reference genomes.')
    parser.add_argument('metadata_file', help='IMG metadata file.')
    parser.add_argument('genome_dir', help='IMG genome directory.')
    parser.add_argument('--num_windows', help='number of windows to sample', type = int, default = 10000)
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)

    args = parser.parse_args()

    deltaTetraDiff = DeltaTetraDiff()
    deltaTetraDiff.run(args.metadata_file, args.genome_dir, args.num_windows, args.threads)
