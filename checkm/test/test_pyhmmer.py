###############################################################################
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

"""Tests of the HMMER functionality provided through pyhmmer."""

import gzip
import os
import shutil
import tempfile
import unittest

import pyhmmer

from checkm.defaultValues import DefaultValues
from checkm.hmmer import HMMERRunner
from checkm.ssuFinder import SSU_Finder

MODEL = 'PF00164.20'  # Ribosom_S12_S23 in the phylogenetic HMMs


def consensus(hmmFile, key=None):
    """Consensus sequence of a HMM (the first HMM if no name or accession is given)."""
    with pyhmmer.plan7.HMMFile(hmmFile) as f:
        for hmm in f:
            if key is None or key in (hmm.name, hmm.accession):
                return hmm.consensus.upper()

    raise KeyError(key)


@unittest.skipUnless(os.path.exists(DefaultValues.PHYLO_HMM_MODELS), 'CheckM reference data is not available')
class VerifyPyhmmer(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm test ')
        self.seqFile = os.path.join(self.tmpDir, 'genes.faa.gz')
        with gzip.open(self.seqFile, 'wt') as f:
            f.write('>gene_1 # a hit\n%s\n' % consensus(DefaultValues.PHYLO_HMM_MODELS, MODEL))
            f.write('>gene_2 # not a hit\n%s\n' % ('MKVLAAGIVGLLLAQWERTYHKLMNPQ' * 4))

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testSearch(self):
        """Verify a search writes the domain table and, if requested, the report with alignments."""
        tableOut = os.path.join(self.tmpDir, 'hits.txt.gz')
        reportOut = os.path.join(self.tmpDir, 'hits.ali.txt.gz')
        HMMERRunner().search(DefaultValues.PHYLO_HMM_MODELS, self.seqFile, tableOut, reportOut, bKeepOutput=True)

        with gzip.open(tableOut, 'rt') as f:
            rows = [line.split() for line in f if line[0] != '#']
        self.assertEqual([(r[0], r[4]) for r in rows], [('gene_1', MODEL)])
        self.assertLess(float(rows[0][6]), 1e-10)

        with gzip.open(reportOut, 'rt') as f:
            report = f.read()
        self.assertIn('Accession:   ' + MODEL, report)
        self.assertIn('>> gene_1  # a hit', report)
        self.assertIn('== domain 1', report)
        self.assertIn('[No hits detected that satisfy reporting thresholds]', report)  # other HMMs
        self.assertTrue(report.endswith('[ok]\n'))

    def testSearchWithoutReport(self):
        """Verify the report is not written unless requested."""
        tableOut = os.path.join(self.tmpDir, 'hits.txt')
        reportOut = os.path.join(self.tmpDir, 'hits.ali.txt')
        HMMERRunner().search(DefaultValues.PHYLO_HMM_MODELS, self.seqFile, tableOut, reportOut, bKeepOutput=False)
        self.assertTrue(os.path.exists(tableOut))
        self.assertFalse(os.path.exists(reportOut))

    def testFetch(self):
        """Verify HMMs are fetched by name or accession, with and without a SSI index, in the order requested."""
        hmmFile = os.path.join(self.tmpDir, 'models.hmm')
        shutil.copyfile(DefaultValues.PHYLO_HMM_MODELS, hmmFile)
        keyFile = os.path.join(self.tmpDir, 'keys.txt')
        with open(keyFile, 'w') as f:
            f.write('TIGR00422\n%s\n' % MODEL)

        HF = HMMERRunner(mode='fetch')
        outputs = []
        for bIndexed in (False, True):
            if bIndexed:
                HF.index(hmmFile)
            self.assertEqual(os.path.exists(hmmFile + '.ssi'), bIndexed)

            out = os.path.join(self.tmpDir, 'fetched_%s.hmm' % bIndexed)
            HF.fetch(hmmFile, keyFile, out, bKeyFile=True)
            with pyhmmer.plan7.HMMFile(out) as f:
                self.assertEqual([hmm.name for hmm in f], ['TIGR00422', 'Ribosom_S12_S23'])
            with open(out, 'rb') as f:
                outputs.append(f.read())

            with self.assertLogs('timestamp', level='ERROR'):
                with self.assertRaises(SystemExit) as cm:
                    HF.fetch(hmmFile, 'NOT_A_MODEL', os.path.join(self.tmpDir, 'missing.hmm'))
            self.assertEqual(cm.exception.code, 1)

        self.assertEqual(outputs[0], outputs[1])

    def testAlign(self):
        """Verify sequences are aligned to a HMM in Pfam format with the mask of consensus columns."""
        modelFile = os.path.join(self.tmpDir, 'model.hmm')
        HMMERRunner(mode='fetch').fetch(DefaultValues.PHYLO_HMM_MODELS, MODEL, modelFile)

        seq = consensus(modelFile)
        unaligned = os.path.join(self.tmpDir, 'unaligned.faa')
        with open(unaligned, 'w') as f:
            f.write('>bin&&gene_1 [e-value=1e-50,score=180.0]\n%s\n' % seq)
            f.write('>bin&&gene_2\n%s\n' % seq[10:-10])

        aligned = os.path.join(self.tmpDir, 'aligned.faa')
        HMMERRunner(mode='align').align(modelFile, unaligned, aligned, outputFormat='Pfam', trim=False)
        with open(aligned) as f:
            lines = f.read().splitlines()

        self.assertEqual(lines[0], '# STOCKHOLM 1.0')
        self.assertIn('#=GS bin&&gene_1 DE [e-value=1e-50,score=180.0]', lines)
        rows = {line.split()[0]: line.split()[1] for line in lines if line and line[0] != '#' and line != '//'}
        mask = [line.split()[2] for line in lines if line.startswith('#=GC RF')][0]
        self.assertEqual(rows['bin&&gene_1'].upper(), seq)
        self.assertEqual(mask, 'x' * len(seq))
        self.assertEqual(len(rows['bin&&gene_2']), len(seq))


@unittest.skipUnless(os.path.exists(os.path.join(DefaultValues.CHECKM_DATA_DIR, 'hmms_ssu', 'SSU_bacteria.hmm')),
                     'CheckM reference data is not available')
class VerifySsuSearch(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm test ')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testNhmmer(self):
        """Verify 16S genes are found on both strands of plain and gzip-compressed sequences."""
        finder = SSU_Finder(threads=2)
        ssu = consensus(finder.bacteriaModelFile).replace('U', 'T')
        flank = 'ACGTTGCAAGCTTGACCATGGTCA' * 20
        revComp = ssu.translate(finder.compl)[::-1]

        for ext, opener in (('fna', open), ('fna.gz', gzip.open)):
            seqFile = os.path.join(self.tmpDir, 'contigs.' + ext)
            with opener(seqFile, 'wt') as f:
                f.write('>forward\n%s\n>reverse\n%s\n>none\n%s\n' % (flank + ssu + flank, flank + revComp, flank * 4))

            prefix = os.path.join(self.tmpDir, 'ssu_' + ext)
            finder._nhmmer(seqFile, 1e-5, prefix, 'bacteria', finder.bacteriaModelFile)
            hits = finder._readHits(prefix + '.bacteria.txt', 'bacteria', 1e-5)

            self.assertEqual(sorted(hits), ['forward', 'reverse'])
            domain, _evalue, start, end, length, revCompStr = hits['forward'][0]
            self.assertEqual((domain, revCompStr), ('bacteria', 'False'))
            self.assertLess(abs(int(start) - (len(flank) + 1)), 40)
            self.assertLess(abs(int(length) - len(ssu)), 80)
            self.assertEqual(hits['reverse'][0][5], 'True')

            with open(prefix + '.bacteria.table.txt') as f:
                self.assertEqual(sorted(line.split()[0] for line in f if line[0] != '#'), ['forward', 'reverse'])
            self.assertEqual([f for f in os.listdir(self.tmpDir) if f.endswith('.tmp')], [])


if __name__ == "__main__":
    unittest.main()
