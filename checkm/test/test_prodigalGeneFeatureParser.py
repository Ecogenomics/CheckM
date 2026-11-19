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

import gzip
import os
import random
import shutil
import tempfile
import unittest

from checkm.prodigal import ProdigalGeneFeatureParser


def writeGFF(path, genes, translationTable=11):
    """Write a Prodigal-style GFF file; genes is a list of (seqId, start, end), 1-based and inclusive."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'wt') as f:
        f.write('##gff-version  3\n')
        for seqId in dict.fromkeys(s for s, _, _ in genes):
            f.write(f'# Sequence Data: seqnum=1;seqlen=100000;seqhdr="{seqId}"\n')
            f.write(f'# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=50.00;transl_table={translationTable};uses_sd=1\n')
            for s, start, end in genes:
                if s == seqId:
                    f.write(f'{seqId}\tProdigal_v2.6.3\tCDS\t{start}\t{end}\t10.0\t+\t0\tID=1_1;partial=00\n')


def maskCodingBases(genes, seqId, start=0, end=None):
    """Reference: count coding bases with a per-base mask (the method used by CheckM <= 1.3.1)."""
    seqGenes = [(s, e) for g, s, e in genes if g == seqId]
    if not seqGenes:
        return 0
    lastCodingBase = max(e for _, e in seqGenes)
    if end is None:
        end = lastCodingBase
    mask = [0] * lastCodingBase
    for s, e in seqGenes:
        for i in range(s - 1, e):
            mask[i] = 1
    return sum(mask[start:end])


class VerifyCodingBases(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def parser(self, genes, name='genes.gff'):
        path = os.path.join(self.tmpDir, name)
        writeGFF(path, genes)
        return ProdigalGeneFeatureParser(path)

    def testKnownCases(self):
        genes = [('a', 1, 10),      # bases 0-9
                 ('a', 5, 20),      # overlaps the first gene: bases 0-19
                 ('a', 21, 30),     # abuts the previous gene: bases 0-29
                 ('a', 41, 50),     # bases 40-49
                 ('a', 43, 45),     # nested in the previous gene
                 ('a', 61, 61),     # a single base: 60
                 ('b', 101, 200)]   # bases 100-199
        p = self.parser(genes)

        self.assertEqual(p.codingBases('a'), 30 + 10 + 1)
        self.assertEqual(p.codingBases('b'), 100)
        self.assertEqual(p.codingBases('missing'), 0)
        self.assertEqual(p.codingBases('a', 0, 10), 10)
        self.assertEqual(p.codingBases('a', 25, 45), 5 + 5)
        self.assertEqual(p.codingBases('a', 30, 40), 0)       # between genes
        self.assertEqual(p.codingBases('a', 44, 46), 2)       # within a gene
        self.assertEqual(p.codingBases('a', 60, 61), 1)
        self.assertEqual(p.codingBases('a', 61, 1000), 0)     # after the last gene
        self.assertEqual(p.codingBases('a', 0, 1000), 41)     # end beyond the last gene
        self.assertEqual(p.codingBases('a', 20, 20), 0)       # empty window
        self.assertEqual(p.codingBases('a', 20, 10), 0)       # end before start
        self.assertEqual(p.codingBases('b', 0, 100), 0)
        self.assertEqual(p.translationTable, 11)

    def testRandomGenesMatchMask(self):
        rng = random.Random(1)
        genes = []
        for seqNum in range(20):
            seqId = f'seq_{seqNum}'
            pos = 1
            for _ in range(rng.randint(1, 60)):
                pos = max(1, pos + rng.randint(-300, 400))     # genes may overlap, nest or abut
                genes.append((seqId, pos, pos + rng.randint(0, 600)))
        rng.shuffle(genes)
        p = self.parser(genes, 'genes.gff.gz')

        seqIds = sorted({g for g, _, _ in genes}) + ['missing']
        for seqId in seqIds:
            self.assertEqual(p.codingBases(seqId), maskCodingBases(genes, seqId), seqId)
            for _ in range(200):
                start = rng.randint(0, 25000)
                end = start + rng.randint(-50, 5000)
                self.assertEqual(p.codingBases(seqId, start, end), maskCodingBases(genes, seqId, start, end),
                                 (seqId, start, end))

    def testReturnsInt(self):
        p = self.parser([('a', 1, 10)])
        self.assertIsInstance(p.codingBases('a'), int)


if __name__ == '__main__':
    unittest.main()
