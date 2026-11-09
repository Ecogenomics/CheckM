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

from checkm.placement import GPlaceRunner, gplaceErrors, parseGPlaceVersion, parseMemorySize


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


class VerifyMemorySize(unittest.TestCase):
    def testParseMemorySize(self):
        """Memory sizes are parsed as by gPlace's --max_memory option (binary units)."""
        GiB = 1024 ** 3
        self.assertEqual(parseMemorySize('64G'), 64 * GiB)
        self.assertEqual(parseMemorySize('64GB'), 64 * GiB)
        self.assertEqual(parseMemorySize('64gib'), 64 * GiB)
        self.assertEqual(parseMemorySize('500M'), 500 * 1024 ** 2)
        self.assertEqual(parseMemorySize('1.5T'), int(1.5 * 1024 ** 4))
        self.assertEqual(parseMemorySize('16'), 16 * GiB)      # a plain number is GiB
        self.assertEqual(parseMemorySize('2048B'), 2048)        # bytes

    def testInvalidMemorySize(self):
        """Invalid memory sizes are rejected."""
        for value in ['', 'G', '64X', '-4G', '64 G B', 'lots', '0']:
            with self.assertRaises(ValueError, msg=value):
                parseMemorySize(value)


class VerifyGPlaceErrors(unittest.TestCase):
    def testGPlaceErrors(self):
        """Error messages are extracted from a gPlace log."""
        tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        try:
            logFile = os.path.join(tmpDir, 'pplacer.out')
            with open(logFile, 'w') as f:
                f.write('[2026-11-09 23:18:01] INFO: Read 128 sequences from concatenated.fasta\n')
                f.write('[2026-11-09 23:18:02] ERROR: --max_memory 4.0 GiB is too small: placing 128 queries '
                        'on this reference package needs --max_memory 5.0G or more.\n')
                f.write('[2026-11-09 23:18:02] ERROR: Controlled exit resulting from an unrecoverable error or warning.\n')
            self.assertEqual(gplaceErrors(logFile),
                             ['--max_memory 4.0 GiB is too small: placing 128 queries on this '
                              'reference package needs --max_memory 5.0G or more.'])

            # a missing log gives no messages
            self.assertEqual(gplaceErrors(os.path.join(tmpDir, 'missing.out')), [])
        finally:
            shutil.rmtree(tmpDir)


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
        self.assertIsNone(runner.maxMemory)

    def testPlaceCmd(self):
        """--max_memory is passed to 'gplace place' only when a memory limit is given."""
        fakeGPlace(self.tmpDir, 'gplace: version 0.1.0 Copyright 2026 Pierre-Alain Chaumeil')
        self._usePath(self.tmpDir + os.pathsep + self.origPath)

        runner = GPlaceRunner(threads=8)
        self.assertEqual(runner._placeCmd('ref.refpkg', 'out.json', 'aln.fasta'),
                         ['gplace', 'place', '-j', '8', '-c', 'ref.refpkg', '-o', 'out.json', 'aln.fasta'])

        runner = GPlaceRunner(threads=8, maxMemory='16G')
        self.assertEqual(runner._placeCmd('ref.refpkg', 'out.json', 'aln.fasta'),
                         ['gplace', 'place', '-j', '8', '-c', 'ref.refpkg', '-o', 'out.json',
                          '--max_memory', '16G', 'aln.fasta'])

    def testInvalidMaxMemory(self):
        """An invalid memory limit is rejected."""
        fakeGPlace(self.tmpDir, 'gplace: version 0.1.0 Copyright 2026 Pierre-Alain Chaumeil')
        self._usePath(self.tmpDir + os.pathsep + self.origPath)
        with self.assertRaises(ValueError):
            GPlaceRunner(threads=1, maxMemory='lots')


if __name__ == "__main__":
    unittest.main()
