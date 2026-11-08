###############################################################################
#
# test_util_coverage.py - unit tests for coverage and profile commands
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
import random
import shutil
import tempfile
import logging
import unittest
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO

try:
    import pysam
except ImportError:
    pysam = None

from checkm.util.seqUtils import writeFasta


def _randomSeq(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def makeBam(bamFile, seqs, readsPerSeq, bNM=True, readLen=50, insert=150):
    """Write a sorted and indexed BAM with perfectly aligned read pairs."""
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'},
              'SQ': [{'SN': seqId, 'LN': len(seq)} for seqId, seq in seqs]}

    rng = random.Random(1)
    unsortedFile = bamFile + '.unsorted.bam'
    with pysam.AlignmentFile(unsortedFile, 'wb', header=header) as out:
        readIndex = 0
        for tid, (_seqId, seq) in enumerate(seqs):
            for _ in range(readsPerSeq):
                start = rng.randint(0, len(seq) - insert)
                mateStart = start + insert - readLen
                for bFirst, pos, matePos in ((True, start, mateStart), (False, mateStart, start)):
                    read = pysam.AlignedSegment()
                    read.query_name = 'r%d' % readIndex
                    read.query_sequence = seq[pos:pos + readLen]
                    read.flag = 0x1 | 0x2 | (0x40 | 0x20 if bFirst else 0x80 | 0x10)
                    read.reference_id = tid
                    read.reference_start = pos
                    read.mapping_quality = 60
                    read.cigartuples = [(0, readLen)]
                    read.next_reference_id = tid
                    read.next_reference_start = matePos
                    read.template_length = insert if bFirst else -insert
                    read.query_qualities = pysam.qualitystring_to_array('I' * readLen)
                    if bNM:
                        read.set_tag('NM', 0)
                    out.write(read)
                readIndex += 1

    pysam.sort('-o', bamFile, unsortedFile)
    pysam.index(bamFile)
    os.remove(unsortedFile)


@unittest.skipIf(pysam is None, 'pysam not installed')
class VerifyCoverage(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

        rng = random.Random(0)
        self.binSeqs = [('c1', _randomSeq(rng, 1000)), ('c2', _randomSeq(rng, 800))]
        self.unbinnedSeqs = [('u%d' % i, _randomSeq(rng, 500)) for i in range(5, 0, -1)]
        self.allSeqs = self.binSeqs + self.unbinnedSeqs

        self.binFile = os.path.join(self.tmpDir, 'bin1.fna')
        writeFasta(dict(self.binSeqs), self.binFile)

        logging.getLogger('timestamp').setLevel(logging.CRITICAL + 1)

    def tearDown(self):
        logging.getLogger('timestamp').setLevel(logging.NOTSET)
        shutil.rmtree(self.tmpDir)

    def _run(self, bamFiles, threads=2):
        from checkm.coverage import Coverage

        outFile = os.path.join(self.tmpDir, 'coverage.tsv')
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            Coverage(threads).run([self.binFile], bamFiles, outFile, False, 0.98, 0.02, 15)

        with open(outFile) as f:
            return [line.rstrip('\n').split('\t') for line in f][1:]

    def testCoverage(self):
        """Coverage is calculated and unbinned sequences are reported in BAM order."""
        bamFile = os.path.join(self.tmpDir, 'sample.bam')
        makeBam(bamFile, self.allSeqs, 10)

        rows = self._run([bamFile])
        self.assertEqual([r[0] for r in rows], ['c1', 'c2', 'u5', 'u4', 'u3', 'u2', 'u1'])
        self.assertEqual([r[1] for r in rows], ['bin1', 'bin1'] + ['unbinned'] * 5)

        for r in rows:
            self.assertEqual(r[3], 'sample')
            self.assertEqual(int(r[5]), 20)  # 10 read pairs per sequence
            self.assertAlmostEqual(float(r[4]), 20 * 50.0 / int(r[2]), places=4)

    def testNoMappedReads(self):
        """A BAM file without reads gives zero coverage instead of an error."""
        bamFile = os.path.join(self.tmpDir, 'empty.bam')
        makeBam(bamFile, self.allSeqs, 0)

        rows = self._run([bamFile])
        self.assertEqual(len(rows), 7)
        for r in rows:
            self.assertEqual(float(r[4]), 0)
            self.assertEqual(int(r[5]), 0)

        # the profile command must handle BAM files without mapped reads
        from checkm.profile import Profile
        profileFile = os.path.join(self.tmpDir, 'profile.tsv')
        with redirect_stdout(StringIO()):
            Profile().run(os.path.join(self.tmpDir, 'coverage.tsv'), profileFile, True)
        with open(profileFile) as f:
            lines = f.readlines()
        self.assertEqual(lines[1].split('\t')[2:4], ['0', '0.0'])

    def testMissingNMTag(self):
        """A BAM file without NM tags is reported as an error."""
        bamFile = os.path.join(self.tmpDir, 'noNM.bam')
        makeBam(bamFile, self.allSeqs, 5, bNM=False)

        with self.assertRaises(SystemExit) as cm:
            self._run([bamFile])
        self.assertEqual(cm.exception.code, 1)

    def testDuplicateBamIds(self):
        """BAM files with the same filename in different directories are reported as an error."""
        bamFiles = []
        for d in ['s1', 's2']:
            os.makedirs(os.path.join(self.tmpDir, d))
            bamFile = os.path.join(self.tmpDir, d, 'aligned.bam')
            makeBam(bamFile, self.allSeqs, 1)
            bamFiles.append(bamFile)

        with self.assertRaises(SystemExit) as cm:
            self._run(bamFiles)
        self.assertEqual(cm.exception.code, 1)

    def testCoverageWindows(self):
        """Coverage of windows (gc_bias_plot) includes a final window ending at the sequence end."""
        from checkm.coverageWindows import CoverageWindows

        bamFile = os.path.join(self.tmpDir, 'sample.bam')
        makeBam(bamFile, self.allSeqs, 10)

        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            profile = CoverageWindows(2).run([self.binFile], bamFile, False, 0.98, 0.02, 100)

        self.assertEqual(sorted(profile.keys()), sorted(seqId for seqId, _ in self.allSeqs))
        self.assertEqual(len(profile['c1'][1]), 10)  # 1000 bp / 100 bp windows
        self.assertEqual(len(profile['c2'][1]), 8)
        self.assertEqual(len(profile['u1'][1]), 5)
        for seqId, seq in self.allSeqs:
            coverage, windowCoverages = profile[seqId]
            self.assertAlmostEqual(coverage, 20 * 50.0 / len(seq))
            self.assertAlmostEqual(sum(windowCoverages) * 100, coverage * len(seq))

    def testCoverageWindowsMissingNMTag(self):
        from checkm.coverageWindows import CoverageWindows

        bamFile = os.path.join(self.tmpDir, 'noNM.bam')
        makeBam(bamFile, self.allSeqs, 5, bNM=False)

        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                CoverageWindows(2).run([self.binFile], bamFile, False, 0.98, 0.02, 100)
        self.assertEqual(cm.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
