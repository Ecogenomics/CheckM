###############################################################################
#
# test_util_genomicSignatures.py - unit tests for tetra and outlier commands
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
import logging
import tempfile
import unittest
from contextlib import redirect_stderr
from io import StringIO

from checkm.genomicSignatures import GenomicSignatures
from checkm.util.seqUtils import writeFasta


class VerifyTetraFile(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        logging.getLogger('timestamp').setLevel(logging.CRITICAL + 1)

    def tearDown(self):
        logging.getLogger('timestamp').setLevel(logging.NOTSET)
        shutil.rmtree(self.tmpDir)

    def testOutputOrder(self):
        """Signatures are written in the order of the input sequences."""
        rng = random.Random(0)
        seqs = {}
        for i in range(200):
            seqs['seq%d' % rng.randint(0, 10 ** 9)] = ''.join(rng.choice('ACGT') for _ in range(rng.randint(10, 3000)))
        seqFile = os.path.join(self.tmpDir, 'seqs.fna')
        writeFasta(seqs, seqFile)

        outFile = os.path.join(self.tmpDir, 'tetra.tsv')
        gs = GenomicSignatures(K=4, threads=4)
        with redirect_stderr(StringIO()):
            gs.calculate(seqFile, outFile)

        with open(outFile) as f:
            seqIds = [line.split('\t')[0] for line in f][1:]
        self.assertEqual(seqIds, list(seqs.keys()))

        sigs = gs.read(outFile)
        for seqId, seq in seqs.items():
            self.assertEqual(list(sigs[seqId]), list(gs.seqSignature(seq)))


if __name__ == "__main__":
    unittest.main()
