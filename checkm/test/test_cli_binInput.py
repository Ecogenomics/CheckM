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

"""Tests for parsing of the bin input (directory or table) used by most commands."""

import os
import shutil
import tempfile
import unittest

from checkm.main import OptionsParser
from checkm.common import checkBinInputExists, checkDirExists, checkFileExists


def writeFile(path, content):
    with open(path, 'w') as f:
        f.write(content)


class TestBinFiles(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        self.parser = OptionsParser()

        self.binDir = os.path.join(self.tmpDir, 'bins')
        os.makedirs(self.binDir)
        self.binA = os.path.join(self.binDir, 'binA.fna')
        self.binB = os.path.join(self.binDir, 'binB.fna')
        writeFile(self.binA, '>s1\nACGT\n')
        writeFile(self.binB, '>s2\nACGT\n')
        self.genesA = os.path.join(self.binDir, 'binA.faa')
        writeFile(self.genesA, '>g1\nMKV\n')
        writeFile(os.path.join(self.binDir, 'empty.fna'), '')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def table(self, content):
        path = os.path.join(self.tmpDir, 'bins.tsv')
        writeFile(path, content)
        return path

    def testDirectory(self):
        binFiles = self.parser.binFiles(self.binDir, 'fna', False)
        self.assertEqual(binFiles, [self.binA, self.binB])

    def testDirectoryNoMatchingExtension(self):
        with self.assertRaises(SystemExit) as cm:
            self.parser.binFiles(self.binDir, 'fasta', False)
        self.assertEqual(cm.exception.code, 1)

    def testMissingInput(self):
        with self.assertRaises(SystemExit) as cm:
            self.parser.binFiles(os.path.join(self.tmpDir, 'missing'), 'fna', False)
        self.assertEqual(cm.exception.code, 1)

    def testTableSkipsBlankAndCommentLines(self):
        path = self.table('# genome_id\tgenome_file\n\nA\t%s\n   \nB\t%s\n\n' % (self.binA, self.binB))
        self.assertEqual(self.parser.binFiles(path, 'fna', False), [self.binA, self.binB])

    def testTableTooFewColumns(self):
        path = self.table('A\t%s\nB\n' % self.binA)
        with self.assertRaises(SystemExit) as cm:
            self.parser.binFiles(path, 'fna', False)
        self.assertEqual(cm.exception.code, 1)

    def testTableCalledGenesNeedsThreeColumns(self):
        path = self.table('A\t%s\n' % self.binA)
        with self.assertRaises(SystemExit) as cm:
            self.parser.binFiles(path, 'fna', True)
        self.assertEqual(cm.exception.code, 1)

        path = self.table('A\t%s\t%s\n' % (self.binA, self.genesA))
        self.assertEqual(self.parser.binFiles(path, 'fna', True), [self.genesA])

    def testTableNoExistingBins(self):
        path = self.table('A\t%s\n' % os.path.join(self.tmpDir, 'missing.fna'))
        with self.assertLogs('timestamp', level='ERROR') as logs:
            with self.assertRaises(SystemExit):
                self.parser.binFiles(path, 'fna', False)
        output = '\n'.join(logs.output)
        self.assertIn('bins input table', output)
        self.assertNotIn('extension (-x)', output)

    def testRedundantBinIdsDifferentExtension(self):
        # x.fna and x.fa both have bin ID 'x' and would overwrite each other's output
        fna = os.path.join(self.tmpDir, 'x.fna')
        fa = os.path.join(self.tmpDir, 'x.fa')
        writeFile(fna, '>s\nACGT\n')
        writeFile(fa, '>s\nACGT\n')
        path = self.table('A\t%s\nB\t%s\n' % (fna, fa))
        with self.assertRaises(SystemExit) as cm:
            self.parser.binFiles(path, 'fna', False)
        self.assertEqual(cm.exception.code, 1)

    def testRedundantBinIdsDifferentDirectory(self):
        otherDir = os.path.join(self.tmpDir, 'other')
        os.makedirs(otherDir)
        otherA = os.path.join(otherDir, 'binA.fna')
        writeFile(otherA, '>s\nACGT\n')
        path = self.table('A\t%s\nB\t%s\n' % (self.binA, otherA))
        with self.assertRaises(SystemExit):
            self.parser.binFiles(path, 'fna', False)

    def testBinsWithGenes(self):
        binFiles = [self.binA, self.binB]
        self.assertEqual(self.parser.binsWithGenes(binFiles, {'binA': {}, 'binB': {}}), binFiles)
        self.assertEqual(self.parser.binsWithGenes(binFiles, {'binB': {}}), [self.binB])
        with self.assertRaises(SystemExit) as cm:
            self.parser.binsWithGenes(binFiles, {})
        self.assertEqual(cm.exception.code, 1)


class TestCheckBinInputExists(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def table(self, content):
        path = os.path.join(self.tmpDir, 'bins.tsv')
        writeFile(path, content)
        return path

    def testDirectory(self):
        checkBinInputExists(self.tmpDir, False)

    def testMissing(self):
        with self.assertRaises(SystemExit):
            checkBinInputExists(os.path.join(self.tmpDir, 'missing'), False)

    def testValidTables(self):
        checkBinInputExists(self.table('\n# comment\nA\ta.fna\n'), False)
        checkBinInputExists(self.table('A\ta.fna\ta.faa\n'), True)

    def testInvalidTables(self):
        with self.assertRaises(SystemExit) as cm:
            checkBinInputExists(self.table('A\n'), False)
        self.assertEqual(cm.exception.code, 1)

        with self.assertRaises(SystemExit):
            checkBinInputExists(self.table('A\ta.fna\n'), True)

        with self.assertRaises(SystemExit):
            checkBinInputExists(self.table('\n'), False)

    def testCheckDirExists(self):
        checkDirExists(self.tmpDir)

        path = self.table('A\ta.fna\n')
        with self.assertRaises(SystemExit):
            checkDirExists(path)

        with self.assertRaises(SystemExit):
            checkDirExists(os.path.join(self.tmpDir, 'missing'))

    def testCheckFileExists(self):
        checkFileExists(self.table('A\ta.fna\n'))
        with self.assertRaises(SystemExit):
            checkFileExists(os.path.join(self.tmpDir, 'missing'))


if __name__ == "__main__":
    unittest.main()
