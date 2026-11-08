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
import tempfile
import unittest

from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.binStatistics import BinStatistics


class VerifyBinStatistics(unittest.TestCase):
    def testGC(self):
        """Verify computation of GC."""
        binStats = BinStatistics(threads=1)

        seqs = {'S1': 'ACgt', 'S2': 'GGgg', 'S3': 'TTtt', 'S4': 'NNNN'}

        seqStats = defaultdict(dict)
        meanGC, _ = binStats.calculateGC(seqs, seqStats)

        self.assertAlmostEqual(seqStats['S1']['GC'], 0.5)
        self.assertAlmostEqual(seqStats['S2']['GC'], 1.0)
        self.assertAlmostEqual(seqStats['S3']['GC'], 0.0)
        self.assertAlmostEqual(seqStats['S4']['GC'], 0.0)

        self.assertAlmostEqual(meanGC, 6.0 / 12.0)

    def testScaffoldLengthStats(self):
        """Verify computation of scaffold length statistics."""
        binStats = BinStatistics(threads=1)

        scaffolds = {'S1': 'ACGT' + DefaultValues.CONTIG_BREAK + 'ACGT', 'S2': 'ACGTACGT', 'S3': 'TTtt'}

        scaffoldStats = defaultdict(dict)
        (maxScaffoldLen, maxContigLen, totalScaffoldBps,
         scaffoldN50, contigN50, scaffoldAvgLen, contigAvgLen,
         numContigs, numAmbiguousBases) = binStats.calculateSeqStats(scaffolds, scaffoldStats)

        self.assertAlmostEqual(scaffoldStats['S1']['Length'], len(DefaultValues.CONTIG_BREAK) + 8)
        self.assertAlmostEqual(scaffoldStats['S1']['Total contig length'], 8)
        self.assertAlmostEqual(scaffoldStats['S1']['# contigs'], 2)

        # scaffold lengths: 18, 8, 4; contig lengths: 4, 4, 8, 4
        self.assertAlmostEqual(maxScaffoldLen, len(DefaultValues.CONTIG_BREAK) + 8)
        self.assertAlmostEqual(maxContigLen, 8)
        self.assertAlmostEqual(totalScaffoldBps, len(DefaultValues.CONTIG_BREAK) + 8 + 8 + 4)
        self.assertAlmostEqual(scaffoldN50, len(DefaultValues.CONTIG_BREAK) + 8)
        self.assertAlmostEqual(contigN50, 4)
        self.assertAlmostEqual(scaffoldAvgLen, 30.0 / 3)
        self.assertAlmostEqual(contigAvgLen, 20.0 / 4)
        self.assertAlmostEqual(numContigs, 4)
        self.assertAlmostEqual(numAmbiguousBases, len(DefaultValues.CONTIG_BREAK))

    def testCodingBases(self):
        """Verify computation of coding bases."""

        binStats = BinStatistics(threads=1)

        aaGenes = {'S1_C1': 'ACGTACGT', 'S1_C2': 'ACGTACGT', 'S3_C1': 'TTtt'}

        codingBasePairs = binStats._calculateCodingBases(aaGenes)

        self.assertAlmostEqual(codingBasePairs, len(aaGenes['S1_C1']) * 3 + len(aaGenes['S1_C2']) * 3 + len(aaGenes['S3_C1']) * 3)

    def _writeBin(self, tmpDir):
        binFile = os.path.join(tmpDir, 'bin1.fna')
        with open(binFile, 'w') as f:
            f.write('>S1\nACGTACGTAC\n>S2\nGGGGCCCC\n')

        return binFile

    def testSequenceStatsWithGenes(self):
        """Verify per-sequence ORF and coding base counts from called genes."""
        binStats = BinStatistics(threads=1)

        with tempfile.TemporaryDirectory() as tmpDir:
            binFile = self._writeBin(tmpDir)

            geneDir = os.path.join(tmpDir, 'bins', 'bin1')
            os.makedirs(geneDir)
            with open(os.path.join(geneDir, DefaultValues.PRODIGAL_AA), 'w') as f:
                f.write('>S1_1\nMK\n>S1_2\nMKV\n')

            seqStats = binStats.sequenceStats(tmpDir, binFile)

        self.assertEqual(seqStats['S1']['# ORFs'], 2)
        self.assertEqual(seqStats['S1']['Coding bases'], (2 + 3) * 3)
        self.assertNotIn('# ORFs', seqStats['S2'])
        self.assertAlmostEqual(seqStats['S2']['GC'], 1.0)
        self.assertEqual(seqStats['S2']['Length'], 8)

    def testSequenceStatsWithoutGenes(self):
        """Verify per-sequence statistics when no called genes are available."""
        binStats = BinStatistics(threads=1)

        with tempfile.TemporaryDirectory() as tmpDir:
            binFile = self._writeBin(tmpDir)
            seqStats = binStats.sequenceStats(tmpDir, binFile)

        for seqId in ['S1', 'S2']:
            self.assertEqual(seqStats[seqId]['# ORFs'], 0)
            self.assertEqual(seqStats[seqId]['Coding bases'], 0)

        self.assertEqual(seqStats['S1']['Length'], 10)
        self.assertAlmostEqual(seqStats['S1']['GC'], 0.5)

if __name__ == "__main__":
    unittest.main()
