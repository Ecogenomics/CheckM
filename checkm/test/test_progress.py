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

import io
import os
import sys
import logging
import tempfile
import shutil
import unittest
import contextlib

from checkm.logger import logger_setup
from checkm.util.progress import progressBar, progressEnabled


class VerifyProgress(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger('timestamp')
        self.origHandlers = list(self.logger.handlers)
        self.origLevel = self.logger.level
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')

    def tearDown(self):
        for handler in list(self.logger.handlers):
            self.logger.removeHandler(handler)
            if handler not in self.origHandlers:
                handler.close()
        for handler in self.origHandlers:
            self.logger.addHandler(handler)
        self.logger.setLevel(self.origLevel)
        shutil.rmtree(self.tmpDir)

    def _setupLogger(self, silent):
        with contextlib.redirect_stdout(io.StringIO()):
            logger_setup(self.tmpDir, 'checkm.log', 'CheckM', 'test', silent)

    def _runBar(self, total, **kwargs):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            progress = progressBar(total, 'Processing bins', 'bin', **kwargs)
            for _ in range(total):
                progress.update(1)
            progress.close()
        return stderr.getvalue()

    def testProgressBar(self):
        """Progress is reported as a tqdm progress bar on stderr."""
        self._setupLogger(silent=False)
        self.assertTrue(progressEnabled())

        output = self._runBar(8)
        lastLine = output.replace('\r', '\n').strip().split('\n')[-1]
        self.assertTrue(lastLine.startswith('    Processing bins: 100%'), lastLine)
        self.assertIn('8/8', lastLine)
        self.assertIn('bin/s', lastLine)

    def testQuiet(self):
        """No progress is reported when console output is suppressed (--quiet)."""
        self._setupLogger(silent=True)
        self.assertFalse(progressEnabled())
        self.assertEqual(self._runBar(8), '')

        # the log file is unaffected
        with open(os.path.join(self.tmpDir, 'checkm.log')) as f:
            self.assertNotIn('Processing bins', f.read())

    def testLoggerLevel(self):
        """No progress is reported if the logger does not report INFO messages."""
        self._setupLogger(silent=False)
        self.logger.setLevel(logging.WARNING)
        self.assertFalse(progressEnabled())
        self.assertEqual(self._runBar(8), '')

    def testDisabled(self):
        """A progress bar can be disabled explicitly."""
        self._setupLogger(silent=False)
        self.assertEqual(self._runBar(8, enabled=False), '')

    def testNoItems(self):
        """A progress bar can be created for no items."""
        self._setupLogger(silent=False)
        self._runBar(0)


if __name__ == "__main__":
    unittest.main()
