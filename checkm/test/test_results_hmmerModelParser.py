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

from checkm.hmmerModelParser import HmmModelParser

MODEL_BODY = """HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.60   4.20   2.80   2.70   3.30   2.90   3.70   2.70   2.60   2.40   3.70   3.10   3.30   3.10   2.90   2.80   2.90   2.60   4.50   3.50
//
"""

# first model has ACC and gathering/trusted/noise thresholds; the second has none of these
HMM_FILE = ("""HMMER3/f [3.1b2 | February 2015]
NAME  modelA
ACC   PF00001.21
DESC  first model
LENG  100
ALPH  amino
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
NSEQ  10
EFFN  1.5
CKSUM 12345
GA    25.00 25.00;
TC    26.00 26.00;
NC    24.00 24.00;
STATS LOCAL MSV       -9.9  0.71
STATS LOCAL VITERBI  -10.7  0.71
STATS LOCAL FORWARD   -4.2  0.71
""" + MODEL_BODY + """HMMER3/f [3.1b2 | February 2015]
NAME  modelB
DESC  second model without accession or thresholds
LENG  50
ALPH  amino
RF    yes
CS    no
MAP   no
NSEQ  5
EFFN  2.5
CKSUM 6789
STATS LOCAL MSV       -8.9  0.71
STATS LOCAL VITERBI  -9.7  0.71
STATS LOCAL FORWARD   -3.2  0.71
""" + MODEL_BODY + """HMMER3/f [3.1b2 | February 2015]
NAME  modelC
ACC   TIGR00001
LENG  75
TC    30.00 30.00;
""" + MODEL_BODY)


class VerifyHmmModelParser(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        self.hmmFile = os.path.join(self.tmpDir, 'test.hmm')
        with open(self.hmmFile, 'w') as f:
            f.write(HMM_FILE)

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testModels(self):
        """Header keys must not leak from one model to the next (simple parse)."""
        models = HmmModelParser(self.hmmFile).models()

        self.assertEqual(sorted(models.keys()), ['PF00001.21', 'TIGR00001', 'modelB'])

        a = models['PF00001.21']
        self.assertEqual(a.name, 'modelA')
        self.assertEqual(a.leng, 100)
        self.assertEqual(a.ga, (25.0, 25.0))
        self.assertEqual(a.tc, (26.0, 26.0))
        self.assertEqual(a.nc, (24.0, 24.0))

        b = models['modelB']
        self.assertEqual(b.name, 'modelB')
        self.assertEqual(b.acc, 'modelB')
        self.assertEqual(b.leng, 50)
        self.assertIsNone(b.ga)
        self.assertIsNone(b.tc)
        self.assertIsNone(b.nc)

        c = models['TIGR00001']
        self.assertEqual(c.name, 'modelC')
        self.assertEqual(c.leng, 75)
        self.assertIsNone(c.ga)
        self.assertEqual(c.tc, (30.0, 30.0))
        self.assertIsNone(c.nc)

    def testParse(self):
        """Header keys must not leak from one model to the next (full parse)."""
        models = list(HmmModelParser(self.hmmFile).parse())
        self.assertEqual([m.acc for m in models], ['PF00001.21', 'modelB', 'TIGR00001'])

        a, b, c = models
        self.assertEqual(a.ga, (25.0, 25.0))
        self.assertEqual(a.nseq, 10)
        self.assertEqual(a.rf, False)
        self.assertEqual(a.map, True)
        self.assertEqual(a.effn, 1.5)
        self.assertEqual(a.stats_local_msv, (-9.9, 0.71))
        self.assertEqual(a.desc, 'first model')

        self.assertIsNone(b.ga)
        self.assertIsNone(b.tc)
        self.assertIsNone(b.nc)
        self.assertEqual(b.rf, True)
        self.assertEqual(b.map, False)
        self.assertEqual(b.nseq, 5)
        self.assertEqual(b.desc, 'second model without accession or thresholds')
        self.assertTrue(b.format.startswith('HMMER3/f'))

        self.assertIsNone(c.ga)
        self.assertEqual(c.tc, (30.0, 30.0))
        self.assertFalse(hasattr(c, 'nseq'))
        self.assertFalse(hasattr(c, 'desc'))


if __name__ == "__main__":
    unittest.main()
