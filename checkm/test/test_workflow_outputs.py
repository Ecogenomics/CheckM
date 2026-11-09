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
from argparse import Namespace

from checkm.main import OptionsParser
from checkm.defaultValues import DefaultValues


class VerifyWorkflowOutputs(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        self.outDir = os.path.join(self.tmpDir, 'output')
        for d in ['bins', 'storage']:
            os.makedirs(os.path.join(self.outDir, d, 'sub'))
        self.genomeDir = os.path.join(self.tmpDir, 'genomes')
        os.makedirs(self.genomeDir)
        self.parser = OptionsParser()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testQualityReportOptions(self):
        """Workflows write the 'qa -o 2 --tab_table' table to <output_dir>/quality_report.tsv."""
        options = Namespace(output_dir=self.outDir)
        self.parser._setQualityReportOptions(options)
        self.assertEqual(options.out_format, 2)
        self.assertTrue(options.bTabTable)
        self.assertEqual(options.file, os.path.join(self.outDir, DefaultValues.QUALITY_REPORT))

    def testKeepTmpFiles(self):
        """The bins and storage directories are kept by default."""
        options = Namespace(output_dir=self.outDir, bin_input=self.genomeDir, bRmTmpFiles=False)
        self.parser._removeTmpFiles(options)
        self.assertTrue(os.path.isdir(os.path.join(self.outDir, 'bins')))
        self.assertTrue(os.path.isdir(os.path.join(self.outDir, 'storage')))

    def testRemoveTmpFiles(self):
        """--rm_tmp_files removes the bins and storage directories only."""
        report = os.path.join(self.outDir, DefaultValues.QUALITY_REPORT)
        open(report, 'w').close()
        options = Namespace(output_dir=self.outDir, bin_input=self.genomeDir, bRmTmpFiles=True)
        self.parser._removeTmpFiles(options)
        self.assertFalse(os.path.exists(os.path.join(self.outDir, 'bins')))
        self.assertFalse(os.path.exists(os.path.join(self.outDir, 'storage')))
        self.assertTrue(os.path.exists(report))

    def testInputGenomesNotRemoved(self):
        """--rm_tmp_files never removes a directory containing the input genomes."""
        binInput = os.path.join(self.outDir, 'bins')
        options = Namespace(output_dir=self.outDir, bin_input=binInput, bRmTmpFiles=True)
        with self.assertLogs('timestamp', level='WARNING'):
            self.parser._removeTmpFiles(options)
        self.assertTrue(os.path.isdir(binInput))
        self.assertFalse(os.path.exists(os.path.join(self.outDir, 'storage')))


if __name__ == "__main__":
    unittest.main()
