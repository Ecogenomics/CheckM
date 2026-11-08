#!/usr/bin/env python

###############################################################################
#
# calculateBoundsDeltaGC.py - find confidence intervals for GC distribution
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
from multiprocessing.connection import wait
import logging

import numpy as np

from checkm.util.seqUtils import readFasta


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


class GenomicSignatures(object):
    def __init__(self, K, threads):
        self.logger = logging.getLogger('timestamp')

        self.K = K
        self.compl = str.maketrans('ACGT', 'TGCA')
        self.kmerCols, self.kmerToCanonicalIndex = self._makeKmerColNames()

        self.totalThreads = threads

    def _makeKmerColNames(self):
        """Work out unique kmers."""

        # determine all mers of a given length
        baseWords = ("A", "C", "G", "T")
        mers = ["A", "C", "G", "T"]
        for _ in range(1, self.K):
            workingList = []
            for mer in mers:
                for char in baseWords:
                    workingList.append(mer + char)
            mers = workingList

        # pare down kmers based on lexicographical ordering
        retList = []
        for mer in mers:
            kmer = self._lexicographicallyLowest(mer)
            if kmer not in retList:
                retList.append(kmer)

        # create mapping from kmers to their canonical order position
        kmerToCanonicalIndex = {}
        for index, kmer in enumerate(retList):
            kmerToCanonicalIndex[kmer] = index
            kmerToCanonicalIndex[self._revComp(kmer)] = index

        return retList, kmerToCanonicalIndex

    def _lexicographicallyLowest(self, seq):
        """Return the lexicographically lowest form of this sequence."""
        rseq = self._revComp(seq)
        if(seq < rseq):
            return seq
        return rseq

    def _revComp(self, seq):
        """Return the reverse complement of a sequence."""
        return seq.translate(self.compl)[::-1]

    def _calculateResults(self, queueIn, queueOut):
        """Calculate genomic signature of sequences in parallel."""
        while True:
            seqId, seq = queueIn.get(block=True, timeout=None)
            if seqId == None:
                break

            sig = self.seqSignature(seq)

            queueOut.put((seqId, sig))

    def _storeResults(self, seqIds, outputFile, writerQueue):
        """Store genomic signatures to file in the order given by seqIds."""

        # write header
        fout = open(outputFile, 'w')
        fout.write('Sequence Id')
        for kmer in self.canonicalKmerOrder():
            fout.write('\t' + kmer)
        fout.write('\n')

        totalSeqs = len(seqIds)
        pendingSigs = {}
        nextSeqIndex = 0
        numProcessedSeq = 0
        while True:
            seqId, sig = writerQueue.get(block=True, timeout=None)
            if seqId == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedSeq += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) sequences.' % (numProcessedSeq, totalSeqs, float(numProcessedSeq) * 100 / totalSeqs)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

            # results arrive in the order workers finish, so
            # write them out in the order of the input sequences
            pendingSigs[seqId] = sig
            while nextSeqIndex < totalSeqs and seqIds[nextSeqIndex] in pendingSigs:
                nextSeqId = seqIds[nextSeqIndex]
                fout.write(nextSeqId)
                fout.write('\t' + '\t'.join(list(map(str, pendingSigs.pop(nextSeqId)))))
                fout.write('\n')
                nextSeqIndex += 1

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

        fout.close()

    def canonicalKmerOrder(self):
        return self.kmerCols

    def seqSignature(self, seq):
        sig = [0] * len(self.kmerCols)

        tmp_seq = seq.upper()

        numMers = len(tmp_seq) - self.K + 1
        for i in range(0, numMers):
            try:
                kmerIndex = self.kmerToCanonicalIndex[tmp_seq[i:i + self.K]]
                sig[kmerIndex] += 1  # Note: a numpy array would be slow here due to this single element increment
            except KeyError:
                # unknown kmer (e.g., contains a N)
                pass

        # normalize; a sequence without any valid kmers (e.g., shorter
        # than K or all N) has a signature of all zeros
        sig = np.array(sig, dtype=float)
        sigSum = np.sum(sig)
        if sigSum > 0:
            sig /= sigSum

        return sig

    def calculate(self, seqFile, outputFile):
        """Calculate genomic signature of each sequence."""

        self.logger.info('Determining tetranucleotide signature of each sequence.')

        # process each sequence in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        seqs = readFasta(seqFile)

        for seqId, seq in seqs.items():
            workerQueue.put((seqId, seq))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        calcProc = [mp.Process(target=self._calculateResults, args=(workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target=self._storeResults, args=(list(seqs.keys()), outputFile, writerQueue))

        try:
            writeProc.start()

            for p in calcProc:
                p.start()

            _waitForWorkers(calcProc, writeProc)

            writerQueue.put((None, None))
            writeProc.join()
        except BaseException:
            # make sure all processes are terminated
            for p in calcProc + [writeProc]:
                if p.is_alive():
                    p.terminate()
            raise

        if any(p.exitcode != 0 for p in calcProc) or writeProc.exitcode != 0:
            self.logger.error('Failed to calculate tetranucleotide signatures for: %s' % seqFile)
            sys.exit(1)

    def distance(self, sig1, sig2):
        return np.sum(np.abs(sig1 - sig2))

    def read(self, tetraProfileFile):
        sig = {}
        with open(tetraProfileFile) as f:
            next(f)
            for line in f:
                lineSplit = line.split('\t')
                sig[lineSplit[0]] = np.array([float(x) for x in lineSplit[1:]])

        return sig
