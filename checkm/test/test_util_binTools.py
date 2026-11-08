###############################################################################
#
# test_util_binTools.py - unit tests for bin statistics edge cases
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
import shutil
import logging
import tempfile
import unittest

import numpy as np

from checkm.binTools import BinTools
from checkm.genomicSignatures import GenomicSignatures
from checkm.unbinned import Unbinned
from checkm.util.seqUtils import writeFasta


class _CodingBases():
    def __init__(self, codingBases):
        self.cb = codingBases

    def codingBases(self, seqId):
        return self.cb[seqId]


class VerifyBinTools(unittest.TestCase):
    def setUp(self):
        self.seqs = {'s1': 'GGCCAATT' * 10, 's2': 'N' * 100, 's3': 'GGGGCCCCAT' * 5}
        self.gs = GenomicSignatures(K=4, threads=1)
        self.tetraSigs = {seqId: self.gs.seqSignature(seq) for seqId, seq in self.seqs.items()}

    def testSeqSignatureNoKmers(self):
        """Sequences without valid tetranucleotides have a zero signature (not NaN)."""
        for seq in ['', 'ACG', 'N' * 50, 'ACGNACGN']:
            sig = self.gs.seqSignature(seq)
            self.assertTrue(np.all(sig == 0))

    def testGcDistAllN(self):
        meanGC, deltaGCs, GCs = BinTools().gcDist(self.seqs)
        self.assertAlmostEqual(meanGC, (40 + 40) / 130.0)
        self.assertTrue(np.isnan(GCs[1]))
        self.assertTrue(np.isnan(deltaGCs[1]))
        self.assertAlmostEqual(GCs[0], 0.5)

        meanGC, deltaGCs, GCs = BinTools().gcDist({'n': 'NNNN'})
        self.assertEqual(meanGC, 0)

    def testCodingDensityEmptySeq(self):
        seqs = {'a': 'A' * 100, 'b': ''}
        meanCD, deltaCDs, CDs = BinTools().codingDensityDist(seqs, _CodingBases({'a': 50, 'b': 0}))
        self.assertAlmostEqual(meanCD, 0.5)
        self.assertAlmostEqual(CDs[0], 0.5)
        self.assertTrue(np.isnan(CDs[1]))

    def testBinTetraSig(self):
        """Sequences without a signature do not make the bin signature NaN."""
        binTools = BinTools()
        binSig = binTools.binTetraSig(self.seqs, self.tetraSigs)
        self.assertTrue(np.all(np.isfinite(binSig)))
        self.assertAlmostEqual(np.sum(binSig), 1.0)

        expected = (self.tetraSigs['s1'] * 80 + self.tetraSigs['s3'] * 50) / 130.0
        self.assertTrue(np.allclose(binSig, expected))

        # NaN signatures from tetra files of older CheckM versions are also ignored
        tetraSigs = dict(self.tetraSigs)
        tetraSigs['s2'] = tetraSigs['s2'] * np.nan
        self.assertTrue(np.allclose(binTools.binTetraSig(self.seqs, tetraSigs), expected))

        meanTD, deltaTDs = binTools.tetraDiffDist(self.seqs, self.gs, self.tetraSigs, binSig)
        self.assertTrue(np.isfinite(meanTD))
        self.assertTrue(np.isnan(deltaTDs[1]))
        self.assertAlmostEqual(meanTD, (deltaTDs[0] + deltaTDs[2]) / 2)


class VerifyUnbinned(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        logging.getLogger('timestamp').setLevel(logging.CRITICAL + 1)

    def tearDown(self):
        logging.getLogger('timestamp').setLevel(logging.NOTSET)
        shutil.rmtree(self.tmpDir)

    def testUnbinned(self):
        binFile = os.path.join(self.tmpDir, 'bin.fna')
        writeFasta({'b1': 'ACGT' * 10}, binFile)
        seqFile = os.path.join(self.tmpDir, 'all.fna')
        writeFasta({'b1': 'ACGT' * 10, 'u1': 'GGGA' * 10, 'u2': 'N' * 40}, seqFile)

        outSeqFile = os.path.join(self.tmpDir, 'unbinned.fna')
        outStatsFile = os.path.join(self.tmpDir, 'unbinned.tsv')
        Unbinned().run([binFile], seqFile, outSeqFile, outStatsFile, 0)

        with open(outStatsFile) as f:
            rows = [line.rstrip('\n').split('\t') for line in f][1:]
        self.assertEqual(rows, [['u1', '40', '75.00'], ['u2', '40', 'NA']])

    def testEmptySeqFile(self):
        binFile = os.path.join(self.tmpDir, 'bin.fna')
        writeFasta({'b1': 'ACGT'}, binFile)
        seqFile = os.path.join(self.tmpDir, 'all.fna')
        open(seqFile, 'w').close()

        Unbinned().run([binFile], seqFile,
                       os.path.join(self.tmpDir, 'unbinned.fna'),
                       os.path.join(self.tmpDir, 'unbinned.tsv'), 0)


if __name__ == "__main__":
    unittest.main()
