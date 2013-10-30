###############################################################################
#
# coverage.py - calculate coverage of all sequences
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
import multiprocessing as mp
import Queue
import time

import pysam

class CoverageStruct():
    def __init__(self, seqLen):
        self.seqLen = seqLen
        self.mappedReads = []
        self.coverage = []

class Coverage():
    """Calculate coverage of all sequences."""
    def __init__(self, threads=1):
        # thready stuff
        self.totalThreads = threads
        self.threadPool = mp.BoundedSemaphore(self.totalThreads)
        self.manager = mp.Manager()
        self.coverageInfo = self.manager.dict()

    def calculate(self, bamFiles, outFile, bPairsOnly=False, bQuiet=False):
        """Calculate coverage of sequences for each BAM file."""

        # process each fasta file
        if not bQuiet:
            print "  Processing %d file(s) with %d threads:" % (len(bamFiles), self.totalThreads)

        fileIndex = 0
        for bamFile in bamFiles:
            if not bQuiet:
                fileIndex += 1
                print '    Processing %s (%d of %d)' % (bamFile, fileIndex, len(bamFiles))

            self.__processBam(bamFile, bPairsOnly, bQuiet)

        if not bQuiet:
            print '  Results written to: %s' % outFile

        fout = open(outFile, 'w')
        for seqId in self.coverageInfo.keys():
            coverageStruct = self.coverageInfo[seqId]
            fout.write(seqId + '\t' + str(coverageStruct.seqLen))
            for i in xrange(len(coverageStruct.coverage)):
                fout.write('\t' + str(coverageStruct.coverage[i]) + '\t' + str(coverageStruct.mappedReads[i]))
            fout.write('\n')
        fout.close()

    def __calculateResults(self, bamIn, queueIn, queueOut, bPairsOnly):
        """Calculate coverage of reference sequences using multiple processes."""
        while True:
          try:
            reference, length = queueIn.get(block=True, timeout=None)
            if reference == None:
                print 'Terminating thread.'
                break

            coverage, mappedReads = self.__processReference(bamIn, reference, length, bPairsOnly)
            queueOut.put((reference, length, coverage, mappedReads))
          except:
            print '[Error] Unknown exception in Coverage.__calculateResults'

    def __storeResults(self, coverageInfo, queue, numReferences, bQuiet):
        """Store coverage information results determined by multiple processes (i.e., to shared memory)."""
        processedRef = 0
        while True:
            try:
                reference, length, coverage, mappedReads = queue.get(block=True, timeout=None)
                if reference == None:
                    break

                if not bQuiet:
                    processedRef += 1
                    processedStr = '      Processed %d of %d (%.2f%%) reference sequences.' % (processedRef, numReferences, float(processedRef)*100/numReferences)
                    sys.stdout.write('%s\r' % processedStr)
                    sys.stdout.flush()

                # this pattern is necessary for the Manager to serialize data
                coverageStruct = coverageInfo.get(reference, CoverageStruct(length))
                coverageStruct.mappedReads.append(mappedReads)
                coverageStruct.coverage.append(coverage)
                coverageInfo[reference] = coverageStruct
            except:
                print '[Error] Unknown exception in Coverage.__storeResults'

        if not bQuiet:
            sys.stdout.write('\n')

    def __processBam(self, bamFile, bPairsOnly, bQuiet):
        """Calculate coverage of sequences in BAM file."""
        bamIn = pysam.Samfile(bamFile, 'rb')

        if not bQuiet:
            print '      Preparing reference sequences'

        references = bamIn.references
        lengths = bamIn.lengths

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for reference, length in zip(references, lengths):
            workerQueue.put([reference, length])

        # Flag used to terminate each worker processes
        # once all sequences have been processes. Note
        # that checking for a Queue.Empty exception
        # doesn't work since the queue may be empty
        # only because data put on the queue is
        # still buffered.
        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        calcProc = [mp.Process(target = self.__calculateResults, args = (bamIn, workerQueue, writerQueue, bPairsOnly)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__storeResults, args = (self.coverageInfo, writerQueue, len(references), bQuiet))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None, None, None))
        writeProc.join()

        bamIn.close()

    def __processReference(self, bamIn, reference, length, bPairsOnly):
        """Calculate coverage of reference sequences."""
        alignedReads = bamIn.fetch(reference, start = 0, end = length, until_eof=True)
        coverage = 0
        mappedReads = 0
        for read in alignedReads:
            if read.is_unmapped or read.is_duplicate or read.is_qcfail or read.is_secondary:
                continue

            if bPairsOnly and not read.is_proper_pair:
                continue

            mappedReads += 1
            coverage += read.rlen

        coverage = float(coverage) / length

        return coverage, mappedReads

