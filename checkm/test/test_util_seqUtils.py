###############################################################################
#
# test_util_seqUtils.py - unit tests for FASTA reading/writing utilities
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
import tempfile
import unittest

from checkm.util.seqUtils import readFasta, isNucleotide, calculateN50


class VerifyReadFasta(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def _write(self, name, content, mode='w'):
        path = os.path.join(self.tmpDir, name)
        with open(path, mode) as f:
            f.write(content)
        return path

    def testNoTrailingNewline(self):
        """Last base is kept when the file does not end with a newline."""
        path = self._write('a.fna', '>s1 desc\nACGT\nACG\n>s2\nGGGTTT')
        seqs = readFasta(path)
        self.assertEqual(seqs['s1'], 'ACGTACG')
        self.assertEqual(seqs['s2'], 'GGGTTT')

    def testWindowsLineEndings(self):
        """CRLF line endings are not included in sequences."""
        path = self._write('a.fna', b'>s1\r\nACGT\r\nAC\r\n', mode='wb')
        self.assertEqual(readFasta(path)['s1'], 'ACGTAC')

    def testHeaderNotTrimmed(self):
        path = self._write('a.fna', '>s1 some desc\nACGT\n')
        self.assertEqual(list(readFasta(path, trimHeader=False).keys()), ['s1 some desc'])


class VerifyIsNucleotide(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def _write(self, content):
        path = os.path.join(self.tmpDir, 'seqs.fa')
        with open(path, 'w') as f:
            f.write(content)
        return path

    def testEmptyRecordSkipped(self):
        """An empty record does not cause a ZeroDivisionError."""
        self.assertTrue(isNucleotide(self._write('>empty\n>s1\nACGTACGTNN\n')))
        self.assertFalse(isNucleotide(self._write('>empty\n>p1\nMKLVWEFPQRS\n')))

    def testOnlyEmptyRecords(self):
        self.assertTrue(isNucleotide(self._write('>e1\n>e2\n')))


class VerifyCalculateN50(unittest.TestCase):
    def testN50(self):
        self.assertEqual(calculateN50([1, 1, 2, 2, 2, 2, 10]), 10)
        self.assertEqual(calculateN50([2, 3, 4, 5, 6]), 5)
        self.assertEqual(calculateN50([7]), 7)

    def testEmpty(self):
        self.assertEqual(calculateN50([]), 0)

    def testInputNotModified(self):
        seqLens = [1, 5, 3, 2]
        calculateN50(seqLens)
        self.assertEqual(seqLens, [1, 5, 3, 2])


if __name__ == "__main__":
    unittest.main()
