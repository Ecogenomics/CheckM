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

"""Tests for the command-line interface (bin/checkm), the help formatter and logging setup."""

import os
import sys
import shutil
import logging
import argparse
import tempfile
import subprocess
import unittest

from checkm.customHelpFormatter import CustomHelpFormatter
from checkm.logger import logger_setup

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CHECKM_SCRIPT = os.path.join(REPO_ROOT, 'bin', 'checkm')


@unittest.skipUnless('CHECKM_DATA_PATH' in os.environ,
                     'CHECKM_DATA_PATH must be set to run bin/checkm without modifying DATA_CONFIG')
@unittest.skipUnless(os.path.exists(CHECKM_SCRIPT), 'bin/checkm not found')
class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def run_checkm(self, *args):
        env = dict(os.environ)
        env['PYTHONPATH'] = REPO_ROOT + os.pathsep + env.get('PYTHONPATH', '')
        return subprocess.run([sys.executable, CHECKM_SCRIPT] + list(args),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, env=env, cwd=self.tmpDir)

    def testHelp(self):
        for flag in ['-h', '--help']:
            p = self.run_checkm(flag)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn('lineage_wf', p.stdout)
            self.assertIn('test         ->', p.stdout)

        p = self.run_checkm()
        self.assertEqual(p.returncode, 0, p.stderr)

    def testBadTmpDir(self):
        p = self.run_checkm('tree', '--tmpdir', os.path.join(self.tmpDir, 'missing'), 'bins', 'out')
        self.assertEqual(p.returncode, 2)
        self.assertIn('must be a valid directory', p.stderr)
        self.assertNotIn('Traceback', p.stderr)

    def testDataHelp(self):
        p = self.run_checkm('data', '-h')
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertNotIn('the the', p.stdout)
        self.assertNotIn('b', [line.strip() for line in p.stdout.splitlines()])

    def testDataBadUsage(self):
        p = self.run_checkm('data', 'setRoot')
        self.assertEqual(p.returncode, 1)

        p = self.run_checkm('data', 'unknownAction')
        self.assertEqual(p.returncode, 1)

    def testMissingBinInput(self):
        p = self.run_checkm('unique', os.path.join(self.tmpDir, 'missing'))
        self.assertEqual(p.returncode, 1)
        self.assertIn('Bin input does not exist', p.stdout)
        self.assertNotIn('Traceback', p.stderr)

    def testUnknownTaxon(self):
        p = self.run_checkm('taxon_set', 'domain', 'NotATaxon', os.path.join(self.tmpDir, 'x.ms'))
        self.assertEqual(p.returncode, 1)

    def testLogMessagesNotDuplicated(self):
        # output directory cannot be created, so no log file can be written
        badOutDir = os.path.join(self.tmpDir, 'file', 'out')
        with open(os.path.join(self.tmpDir, 'file'), 'w') as f:
            f.write('not a directory')

        p = self.run_checkm('tree', os.path.join(self.tmpDir, 'missing'), badOutDir)
        self.assertEqual(p.returncode, 1)
        self.assertEqual(p.stdout.count('Bin input does not exist'), 1)


class TestCustomHelpFormatter(unittest.TestCase):
    def parser(self):
        return argparse.ArgumentParser('prog', formatter_class=CustomHelpFormatter)

    def testDefaultAdded(self):
        p = self.parser()
        p.add_argument('-t', '--threads', type=int, default=4, help='number of threads')
        self.assertIn('number of threads (default: 4)', p.format_help())

    def testExplicitDefaultPlaceholder(self):
        p = self.parser()
        p.add_argument('-t', '--threads', type=int, default=4, help='threads [%(default)s]')
        self.assertIn('threads [4]', p.format_help())
        self.assertNotIn('(default: 4)', p.format_help())

    def testHelpString(self):
        p = self.parser()
        action = p.add_argument('-x', default='fna', help='extension %(default)s')
        formatter = CustomHelpFormatter('prog')
        self.assertEqual(formatter._get_help_string(action), 'extension %(default)s')


class TestLoggerSetup(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger('timestamp')
        self.savedHandlers = list(self.logger.handlers)
        self.savedLevel = self.logger.level
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        for h in list(self.logger.handlers):
            self.logger.removeHandler(h)
            h.close()
        for h in self.savedHandlers:
            self.logger.addHandler(h)
        self.logger.setLevel(self.savedLevel)
        shutil.rmtree(self.tmpDir)

    def testRepeatedSetupDoesNotDuplicateHandlers(self):
        logger_setup(self.tmpDir, 'checkm.log', 'CheckM', '0.0.0', True)
        logger_setup(self.tmpDir, 'checkm.log', 'CheckM', '0.0.0', True)
        self.assertEqual(len(self.logger.handlers), 2)  # stdout + log file

    def testUnusableLogDir(self):
        notADir = os.path.join(self.tmpDir, 'file')
        with open(notADir, 'w') as f:
            f.write('x')

        logger_setup(os.path.join(notADir, 'out'), 'checkm.log', 'CheckM', '0.0.0', True)
        self.assertEqual(len(self.logger.handlers), 1)  # stdout only


if __name__ == "__main__":
    unittest.main()
