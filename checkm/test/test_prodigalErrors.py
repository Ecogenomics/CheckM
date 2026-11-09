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

from checkm.prodigal import ProdigalRunner


FAKE_PRODIGAL = """#!/bin/sh
# fake Prodigal: succeeds for -h, otherwise reports an error and fails
for arg in "$@"; do
    if [ "$arg" = "-h" ]; then
        exit 0
    fi
done
echo "" >&2
echo "Error:  simulated failure" >&2
echo "  on two lines." >&2
exit 18
"""


class VerifyProdigalErrors(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')

        # fake prodigal placed first on the path
        binDir = os.path.join(self.tmpDir, 'bin')
        os.makedirs(binDir)
        fakeProdigal = os.path.join(binDir, 'prodigal')
        with open(fakeProdigal, 'w') as f:
            f.write(FAKE_PRODIGAL)
        os.chmod(fakeProdigal, os.stat(fakeProdigal).st_mode | stat.S_IXUSR)

        self.origPath = os.environ['PATH']
        os.environ['PATH'] = binDir + os.pathsep + self.origPath

        # genome large enough for '-p single' to be tried before '-p meta'
        self.genomeFile = os.path.join(self.tmpDir, 'genome.fna')
        with open(self.genomeFile, 'w') as f:
            f.write('>contig_1\n' + 'ACGTTGCAAG' * 12000 + '\n')

        self.outDir = os.path.join(self.tmpDir, 'out')
        os.makedirs(self.outDir)

    def tearDown(self):
        os.environ['PATH'] = self.origPath
        shutil.rmtree(self.tmpDir)

    def testErrorOutputReported(self):
        """Prodigal's error output is included in the retry warning and in the error message."""
        prodigal = ProdigalRunner(self.outDir)

        with self.assertLogs('timestamp', level='WARNING') as cm:
            with self.assertRaises(SystemExit) as exitCm:
                prodigal.run(self.genomeFile, bNucORFs=False)

        self.assertEqual(exitCm.exception.code, 1)

        warnings = [r.getMessage() for r in cm.records if r.levelname == 'WARNING']
        errors = [r.getMessage() for r in cm.records if r.levelname == 'ERROR']

        self.assertTrue(any("'-p single' failed (exit code 18: Error: simulated failure on two lines.)" in w
                            for w in warnings), warnings)
        self.assertEqual(len(errors), 1)
        self.assertIn("Prodigal failed with '-p meta' (exit code 18: Error: simulated failure on two lines.)", errors[0])
        self.assertIn(self.genomeFile, errors[0])

    def testErrorOutputFormatting(self):
        """Error output is reported on a single line and truncated if very long."""
        self.assertEqual(ProdigalRunner._errorOutput(b'\nError:  no input\n sequences.\n'), 'Error: no input sequences.')
        self.assertEqual(ProdigalRunner._errorOutput(b''), '')
        self.assertEqual(ProdigalRunner._errorOutput(b'\xff bad bytes'), '� bad bytes')

        longOutput = ProdigalRunner._errorOutput(b'x' * 5000)
        self.assertEqual(len(longOutput), 1000 + len(' ...'))
        self.assertTrue(longOutput.endswith(' ...'))

    def testFailureStr(self):
        """Failures are described by the return code and any error output."""
        self.assertEqual(ProdigalRunner._failureStr(18, 'Error: no input sequences to analyze.'),
                         'exit code 18: Error: no input sequences to analyze.')
        self.assertEqual(ProdigalRunner._failureStr(5, ''), 'exit code 5')
        self.assertEqual(ProdigalRunner._failureStr(-9, ''), 'signal SIGKILL')


if __name__ == "__main__":
    unittest.main()
