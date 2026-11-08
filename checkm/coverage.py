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
import os
import multiprocessing as mp
from multiprocessing.connection import wait
import logging
import ntpath
from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.common import reassignStdOut, restoreStdOut, binIdFromFilename
from checkm.util.seqUtils import readFasta

from numpy import mean, sqrt


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


class CoverageStruct():
    def __init__(self, seqLen, mappedReads, coverage):
        self.seqLen = seqLen
        self.mappedReads = mappedReads
        self.coverage = coverage


class Coverage():
    """Calculate coverage of all sequences."""

    def __init__(self, threads):
        """Initialization."""

        # pysam can be a challenge to install and most users
        # will not use this functionality so importing pysam here
        import pysam
        self.pysam = pysam

        self.logger = logging.getLogger('timestamp')
        self.totalThreads = threads

    def run(self, binFiles, bamFiles, outFile, bAllReads, minAlignPer, maxEditDistPer, minQC):
        """Calculate coverage of sequences for each BAM file."""

        # determine bin assignment of each sequence
        self.logger.info('Determining bin assignment of each sequence.')

        # sequences are tracked as (bin ID, sequence ID) pairs since bins
        # (e.g., from different assemblies) may contain sequences with the
        # same ID, which previously caused bins to be dropped (issue #386)
        binSeqs = []
        seqIdToBinIds = defaultdict(list)
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)

            seqs = readFasta(binFile)
            for seqId, seq in seqs.items():
                binSeqs.append((binId, seqId, len(seq)))
                seqIdToBinIds[seqId].append(binId)

        dupSeqIds = [seqId for seqId, binIds in seqIdToBinIds.items() if len(binIds) > 1]
        if dupSeqIds:
            self.logger.warning('Identified %d sequence IDs present in multiple bins (e.g., %s in bins %s).'
                                % (len(dupSeqIds), dupSeqIds[0], ', '.join(seqIdToBinIds[dupSeqIds[0]])))
            self.logger.warning('BAM files identify sequences only by ID. Mapped reads are assigned to a duplicated '
                                'sequence only if its length matches the BAM reference sequence, and will be assigned '
                                'to all bins containing a sequence with the same ID and length.')

        # process each BAM file
        self.logger.info("Processing %d file(s) with %d threads.\n" % (len(bamFiles), self.totalThreads))

        # BAM files are identified in the output by their filename without
        # path or extension so these must be unique
        bamIdToFiles = defaultdict(list)
        for bamFile in bamFiles:
            bamIdToFiles[binIdFromFilename(bamFile)].append(bamFile)

        dupBamIds = [bamId for bamId, files in bamIdToFiles.items() if len(files) > 1]
        if dupBamIds:
            self.logger.error('BAM files must have unique filenames (ignoring path and extension) as these are used to identify them in the output.')
            for bamId in dupBamIds:
                self.logger.error('BAM ID %s is shared by: %s' % (bamId, ', '.join(bamIdToFiles[bamId])))
            sys.exit(1)

        # make sure all BAM files are sorted
        self.numFiles = len(bamFiles)
        for bamFile in bamFiles:
            if not os.path.exists(bamFile + '.bai'):
                self.logger.error('BAM file is either unsorted or not indexed: ' + bamFile + '\n')
                sys.exit(1)

        # calculate coverage of each BAM file
        coverageInfo = {}
        refSeqOrder = []
        numFilesStarted = 0
        for bamFile in bamFiles:
            numFilesStarted += 1
            self.logger.info('Processing %s (%d of %d):' % (ntpath.basename(bamFile), numFilesStarted, len(bamFiles)))

            with self.pysam.AlignmentFile(bamFile, 'rb') as bam:
                refSeqOrder.extend(bam.references)

            with mp.Manager() as manager:
                coverageInfo[bamFile] = dict(self._processBam(bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, manager.dict()))

        # redirect output
        self.logger.info('Writing coverage information to file.')
        oldStdOut = reassignStdOut(outFile)

        header = 'Sequence Id\tBin Id\tSequence length (bp)'
        for bamFile in bamFiles:
            header += '\tBam Id\tCoverage\tMapped reads'

        print(header)

        # write coverage stats for all binned sequences to file
        for binId, seqId, seqLen in binSeqs:
            bDupSeqId = len(seqIdToBinIds[seqId]) > 1

            rowStr = seqId + '\t' + binId + '\t' + str(seqLen)
            for bamFile in bamFiles:
                bamId = binIdFromFilename(bamFile)

                covInfo = coverageInfo[bamFile].get(seqId)
                if covInfo is not None and bDupSeqId and covInfo.seqLen != seqLen:
                    # BAM reference is a different sequence with the same ID
                    covInfo = None

                if covInfo is not None:
                    rowStr += '\t%s\t%f\t%d' % (bamId, covInfo.coverage, covInfo.mappedReads)
                else:
                    rowStr += '\t%s\t%f\t%d' % (bamId, 0, 0)

            print(rowStr)

        # write coverage stats for unbinned sequences in the BAM files to file,
        # with sequences in the order they appear in the BAM headers
        unbinnedSeqLens = {}
        for seqId in refSeqOrder:
            if seqId in seqIdToBinIds or seqId in unbinnedSeqLens:
                continue

            for bamFile in bamFiles:
                covInfo = coverageInfo[bamFile].get(seqId)
                if covInfo is not None:
                    unbinnedSeqLens[seqId] = covInfo.seqLen
                    break

        for seqId, seqLen in unbinnedSeqLens.items():
            rowStr = seqId + '\t' + DefaultValues.UNBINNED + '\t' + str(seqLen)
            for bamFile in bamFiles:
                bamId = binIdFromFilename(bamFile)

                covInfo = coverageInfo[bamFile].get(seqId)
                if covInfo is not None:
                    rowStr += '\t%s\t%f\t%d' % (bamId, covInfo.coverage, covInfo.mappedReads)
                else:
                    rowStr += '\t%s\t%f\t%d' % (bamId, 0, 0)

            print(rowStr)

        # restore stdout
        restoreStdOut(outFile, oldStdOut)

    def _processBam(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, coverageInfo):
        """Calculate coverage of sequences in BAM file."""

        # determine coverage for each reference sequence
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        bamfile = self.pysam.Samfile(bamFile, 'rb')
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

        workerProc = [mp.Process(target=self.__workerThread, args=(bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target=self.__writerThread, args=(coverageInfo, len(refSeqIds), writerQueue))

        try:
            writeProc.start()

            for p in workerProc:
                p.start()

            _waitForWorkers(workerProc, writeProc)

            writerQueue.put((None, None, None, None, None, None, None, None, None, None, None))
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

    def __workerThread(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            seqIds, seqLens = queueIn.get(block=True, timeout=None)
            if seqIds == None:
                break

            bamfile = self.pysam.Samfile(bamFile, 'rb')

            for seqId, seqLen in zip(seqIds, seqLens):
                numReads = 0
                numMappedReads = 0
                numDuplicates = 0
                numSecondary = 0
                numFailedQC = 0
                numFailedAlignLen = 0
                numFailedEditDist = 0
                numFailedProperPair = 0
                coverage = 0

                for read in bamfile.fetch(seqId, 0, seqLen):
                    numReads += 1

                    if read.is_unmapped:
                        pass
                    elif read.is_duplicate:
                        numDuplicates += 1
                    elif read.is_secondary or read.is_supplementary:
                        numSecondary += 1
                    elif read.is_qcfail or read.mapping_quality < minQC:
                        numFailedQC += 1
                    elif read.query_alignment_length < minAlignPer * read.query_length:
                        numFailedAlignLen += 1
                    elif not read.has_tag('NM'):
                        self.logger.error('Read %s in BAM file %s has no NM (edit distance) tag. '
                                          'CheckM requires BAM files with NM tags (e.g., from BWA or Bowtie2, '
                                          'or added with samtools calmd).' % (read.query_name, bamFile))
                        sys.exit(1)
                    elif read.get_tag('NM') > maxEditDistPer * read.query_length:
                        numFailedEditDist += 1
                    elif not bAllReads and not read.is_proper_pair:
                        numFailedProperPair += 1
                    else:
                        numMappedReads += 1

                        # Note: the alignment length (query_alignment_length) is used instead of the
                        # read length (query_length) as this brings the calculated coverage
                        # in line with 'samtools depth' (at least when the min
                        # alignment length and edit distance thresholds are zero).
                        coverage += read.query_alignment_length

                if seqLen > 0:
                    coverage = float(coverage) / seqLen
                else:
                    coverage = 0.0

                queueOut.put((seqId, seqLen, coverage, numReads,
                                numDuplicates, numSecondary, numFailedQC,
                                numFailedAlignLen, numFailedEditDist,
                                numFailedProperPair, numMappedReads))

            bamfile.close()

    def __writerThread(self, coverageInfo, numRefSeqs, writerQueue):
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
            seqId, seqLen, coverage, numReads, numDuplicates, numSecondary, numFailedQC, numFailedAlignLen, numFailedEditDist, numFailedProperPair, numMappedReads = writerQueue.get(block=True, timeout=None)
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

            coverageInfo[seqId] = CoverageStruct(seqLen=seqLen, mappedReads=numMappedReads, coverage=coverage)

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

    def binProfiles(self, coverageFile):
        """Read coverage information for each bin."""
        binCoverages = defaultdict(lambda: defaultdict(list))
        binStats = defaultdict(dict)

        bHeader = True
        for line in open(coverageFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            binId = lineSplit[1]
            seqLen = int(lineSplit[2])

            # calculate mean coverage (weighted by scaffold length)
            # for each bin under each BAM file
            for i in range(3, len(lineSplit), 3):
                bamId = lineSplit[i]
                coverage = float(lineSplit[i + 1])
                binCoverages[binId][bamId].append(coverage)

                if bamId not in binStats[binId]:
                    binStats[binId][bamId] = [0, 0]

                binLength = binStats[binId][bamId][0] + seqLen
                weight = float(seqLen) / binLength
                meanBinCoverage = coverage * weight + binStats[binId][bamId][1] * (1 - weight)

                binStats[binId][bamId] = [binLength, meanBinCoverage]

        profiles = defaultdict(dict)
        for binId in binStats:
            for bamId, stats in binStats[binId].items():
                binLength, meanBinCoverage = stats
                coverages = binCoverages[binId][bamId]

                varCoverage = 0
                if len(coverages) > 1:
                    varCoverage = mean(list(map(lambda x: (x - meanBinCoverage) ** 2, coverages)))

                profiles[binId][bamId] = [meanBinCoverage, sqrt(varCoverage)]

        return profiles
