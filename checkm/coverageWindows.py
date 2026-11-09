###############################################################################
#
# coverageWindows.py - calculate coverage of windows within sequences
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
import multiprocessing as mp
from multiprocessing.connection import wait
import logging

from numpy import zeros


def _waitForWorkers(workerProcs, writerProc):
    """Wait for worker processes to finish.

    Workers are terminated if the writer process exits
    early, since workers can block on exit until the
    data they put on the writer queue has been consumed.
    """
    pending = list(workerProcs)
    while pending:
        wait([p.sentinel for p in pending] + [writerProc.sentinel])
        pending = [p for p in pending if p.is_alive()]
        if pending and not writerProc.is_alive():
            for p in pending:
                p.terminate()
            for p in pending:
                p.join()
            break

    for p in workerProcs:
        p.join()


class MissingNMTagError(Exception):
    """Read is missing the NM (edit distance) tag."""
    pass


class ReadLoader:
    """Count aligned reads returned by pysam's fetch and calculate per-base coverage."""

    def __init__(self, refLength, bAllReads, minAlignPer, maxEditDistPer):
        """Initialization."""

        self.bAllReads = bAllReads
        self.minAlignPer = minAlignPer
        self.maxEditDistPer = maxEditDistPer

        self.numReads = 0
        self.numMappedReads = 0
        self.numDuplicates = 0
        self.numSecondary = 0
        self.numFailedQC = 0
        self.numFailedAlignLen = 0
        self.numFailedEditDist = 0
        self.numFailedProperPair = 0

        self.coverage = zeros(refLength)

    def __call__(self, read):
        self.numReads += 1

        if read.is_unmapped:
            pass
        elif read.is_duplicate:
            self.numDuplicates += 1
        elif read.is_secondary or read.is_supplementary:
            self.numSecondary += 1
        elif read.is_qcfail:
            self.numFailedQC += 1
        elif read.query_alignment_length < self.minAlignPer * read.query_length:
            self.numFailedAlignLen += 1
        elif not read.has_tag('NM'):
            raise MissingNMTagError(read.query_name)
        elif read.get_tag('NM') > self.maxEditDistPer * read.query_length:
            self.numFailedEditDist += 1
        elif not self.bAllReads and not read.is_proper_pair:
            self.numFailedProperPair += 1
        else:
            self.numMappedReads += 1

            # Note: the aligned reference bases are counted instead of the
            # read length (query_length) as this brings the calculated coverage
            # in line with 'samtools depth' (at least when the min
            # alignment length and edit distance thresholds are zero).
            self.coverage[read.reference_start:read.reference_end] += 1.0


class CoverageWindows():
    """Calculate coverage of all sequences."""
    def __init__(self, threads):
        """Initialization."""

        # pysam can be a challenge to install and most users
        # will not use this functionality so importing pysam here
        import pysam
        self.pysam = pysam

        self.logger = logging.getLogger('timestamp')

        self.totalThreads = threads

    def run(self, binFiles, bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize):
        """Calculate coverage of full sequences and windows."""

        # make sure BAM file is sorted
        if not os.path.exists(bamFile + '.bai'):
            self.logger.error('BAM file is not sorted: ' + bamFile + '\n')
            sys.exit(1)

        # calculate coverage of each BAM file
        self.logger.info('Calculating coverage of windows.')
        with mp.Manager() as manager:
            coverageInfo = dict(self._processBam(bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, manager.dict()))

        return coverageInfo

    def _processBam(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, coverageInfo):
        """Calculate coverage of sequences in BAM file."""

        # determine coverage for each reference sequence
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        bamfile = self.pysam.AlignmentFile(bamFile, 'rb')
        refSeqIds = bamfile.references
        refSeqLens = bamfile.lengths

        # populate each thread with reference sequence to process
        # Note: reference sequences are sorted by number of mapped reads
        # so it is important to distribute reads in a sensible way to each
        # of the threads
        refSeqLists = [[] for _ in range(self.totalThreads)]
        refLenLists = [[] for _ in range(self.totalThreads)]

        threadIndex = 0
        incDir = 1
        for refSeqId, refLen in zip(refSeqIds, refSeqLens):
            refSeqLists[threadIndex].append(refSeqId)
            refLenLists[threadIndex].append(refLen)

            threadIndex += incDir
            if threadIndex == self.totalThreads:
                threadIndex = self.totalThreads - 1
                incDir = -1
            elif threadIndex == -1:
                threadIndex = 0
                incDir = 1

        for i in range(self.totalThreads):
            workerQueue.put((refSeqLists[i], refLenLists[i]))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        bamfile.close()

        workerProc = [mp.Process(target=self._workerThread, args=(bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target=self._writerThread, args=(coverageInfo, len(refSeqIds), writerQueue))

        try:
            writeProc.start()

            for p in workerProc:
                p.start()

            _waitForWorkers(workerProc, writeProc)

            writerQueue.put((None, None, None, None, None, None, None, None, None, None, None, None))
            writeProc.join()
        except BaseException:
            # make sure all processes are terminated
            for p in workerProc + [writeProc]:
                if p.is_alive():
                    p.terminate()
            raise

        if any(p.exitcode != 0 for p in workerProc) or writeProc.exitcode != 0:
            self.logger.error('Failed to calculate coverage for BAM file: %s' % bamFile)
            sys.exit(1)

        return coverageInfo

    def _workerThread(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            seqIds, seqLens = queueIn.get(block=True, timeout=None)
            if seqIds == None:
                break

            bamfile = self.pysam.AlignmentFile(bamFile, 'rb')

            for seqId, seqLen in zip(seqIds, seqLens):
                readLoader = ReadLoader(seqLen, bAllReads, minAlignPer, maxEditDistPer)
                try:
                    for read in bamfile.fetch(seqId, 0, seqLen):
                        readLoader(read)
                except MissingNMTagError as e:
                    self.logger.error('Read %s in BAM file %s has no NM (edit distance) tag. '
                                      'CheckM requires BAM files with NM tags (e.g., from BWA or Bowtie2, '
                                      'or added with samtools calmd).' % (str(e), bamFile))
                    sys.exit(1)

                # coverage of each complete window; a partial window
                # at the end of the sequence is ignored
                start = 0
                end = windowSize
                windowCoverages = []
                while end <= seqLen:
                    windowCoverages.append(sum(readLoader.coverage[start:end]) / windowSize)

                    start = end
                    end += windowSize

                if seqLen > 0:
                    coverage = float(sum(readLoader.coverage)) / seqLen
                else:
                    coverage = 0.0

                queueOut.put((seqId, seqLen, coverage, windowCoverages, readLoader.numReads,
                                readLoader.numDuplicates, readLoader.numSecondary, readLoader.numFailedQC,
                                readLoader.numFailedAlignLen, readLoader.numFailedEditDist,
                                readLoader.numFailedProperPair, readLoader.numMappedReads))

            bamfile.close()

    def _writerThread(self, coverageInfo, numRefSeqs, writerQueue):
        """Store or write results of worker threads in a single thread."""
        totalReads = 0
        totalDuplicates = 0
        totalSecondary = 0
        totalFailedQC = 0
        totalFailedAlignLen = 0
        totalFailedEditDist = 0
        totalFailedProperPair = 0
        totalMappedReads = 0

        processedRefSeqs = 0
        while True:
            seqId, seqLen, coverage, windowCoverages, numReads, numDuplicates, numSecondary, numFailedQC, numFailedAlignLen, numFailedEditDist, numFailedProperPair, numMappedReads = writerQueue.get(block=True, timeout=None)
            if seqId == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedRefSeqs += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) reference sequences.' % (processedRefSeqs, numRefSeqs, float(processedRefSeqs) * 100 / numRefSeqs)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

                totalReads += numReads
                totalDuplicates += numDuplicates
                totalSecondary += numSecondary
                totalFailedQC += numFailedQC
                totalFailedAlignLen += numFailedAlignLen
                totalFailedEditDist += numFailedEditDist
                totalFailedProperPair += numFailedProperPair
                totalMappedReads += numMappedReads

            coverageInfo[seqId] = [coverage, windowCoverages]

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

            def perc(count):
                return float(count) * 100 / totalReads if totalReads > 0 else 0.0

            print ('')
            print ('    # total reads: %d' % totalReads)
            print ('      # properly mapped reads: %d (%.1f%%)' % (totalMappedReads, perc(totalMappedReads)))
            print ('      # duplicate reads: %d (%.1f%%)' % (totalDuplicates, perc(totalDuplicates)))
            print ('      # secondary reads: %d (%.1f%%)' % (totalSecondary, perc(totalSecondary)))
            print ('      # reads failing QC: %d (%.1f%%)' % (totalFailedQC, perc(totalFailedQC)))
            print ('      # reads failing alignment length: %d (%.1f%%)' % (totalFailedAlignLen, perc(totalFailedAlignLen)))
            print ('      # reads failing edit distance: %d (%.1f%%)' % (totalFailedEditDist, perc(totalFailedEditDist)))
            print ('      # reads not properly paired: %d (%.1f%%)' % (totalFailedProperPair, perc(totalFailedProperPair)))
            print ('')
