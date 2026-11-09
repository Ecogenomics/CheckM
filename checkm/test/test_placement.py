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
import stat
import tempfile
import unittest

from checkm.placement import GPlaceRunner, parseGPlaceVersion


def fakeGPlace(binDir, versionOutput, exitCode=0):
    """Create a fake gplace command that reports the given version."""
    path = os.path.join(binDir, 'gplace')
    with open(path, 'w') as f:
        f.write('#!/bin/sh\n')
        f.write("echo '%s'\n" % versionOutput)
        f.write('exit %d\n' % exitCode)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


class VerifyGPlaceVersion(unittest.TestCase):
    def testParseVersion(self):
        """Versions are parsed from the output of 'gplace --version'."""
        self.assertEqual(parseGPlaceVersion('gplace: version 0.1.0 Copyright 2026 Pierre-Alain Chaumeil'), (0, 1, 0))
        self.assertEqual(parseGPlaceVersion('gplace: version 0.10.2'), (0, 10, 2))
        self.assertEqual(parseGPlaceVersion('gplace: version 1.2 Copyright'), (1, 2))
        self.assertIsNone(parseGPlaceVersion('usage: gplace ...'))
        self.assertIsNone(parseGPlaceVersion(''))

        # versions are compared numerically
        self.assertTrue(parseGPlaceVersion('gplace: version 0.10.0') > (0, 9, 0))


class VerifyGPlaceCheck(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        self.origPath = os.environ['PATH']

    def tearDown(self):
        os.environ['PATH'] = self.origPath
        shutil.rmtree(self.tmpDir)

    def _usePath(self, path):
        os.environ['PATH'] = path

    def testMissing(self):
        """CheckM exits with an error if gplace is not on the system path."""
        self._usePath(self.tmpDir)  # an empty directory
        with self.assertLogs('timestamp', level='ERROR') as cm:
            with self.assertRaises(SystemExit) as exitCm:
                GPlaceRunner(threads=1)
        self.assertEqual(exitCm.exception.code, 1)
        self.assertTrue(any('pip install gtdb-gplace' in msg for msg in cm.output))

    def testTooOld(self):
        """CheckM exits with an error if gplace is older than 0.1.0."""
        fakeGPlace(self.tmpDir, 'gplace: version 0.0.9 Copyright 2026')
        self._usePath(self.tmpDir + os.pathsep + self.origPath)
        with self.assertLogs('timestamp', level='ERROR') as cm:
            with self.assertRaises(SystemExit):
                GPlaceRunner(threads=1)
        self.assertTrue(any('gPlace 0.0.9 was found' in msg for msg in cm.output))

    def testUnknownVersion(self):
        """CheckM exits with an error if the version of gplace can't be determined."""
        fakeGPlace(self.tmpDir, 'something unexpected', exitCode=2)
        self._usePath(self.tmpDir + os.pathsep + self.origPath)
        with self.assertLogs('timestamp', level='ERROR'):
            with self.assertRaises(SystemExit):
                GPlaceRunner(threads=1)

    def testSupportedVersion(self):
        """A supported version of gplace is accepted."""
        fakeGPlace(self.tmpDir, 'gplace: version 0.1.0 Copyright 2026 Pierre-Alain Chaumeil')
        self._usePath(self.tmpDir + os.pathsep + self.origPath)
        runner = GPlaceRunner(threads=96)
        self.assertEqual(runner.version, (0, 1, 0))
        self.assertEqual(runner.numThreads, 96)


if __name__ == "__main__":
    unittest.main()
