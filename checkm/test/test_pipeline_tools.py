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

"""Tests for running external tools, worker failures and bins without genes.

Tests that require external programs (prodigal, HMMER) or CheckM reference
data are skipped if these are not available.
"""

import os
import sys
import shutil
import logging
import tempfile
import subprocess
import unittest
from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.hmmer import HMMERRunner, returnCodeStr, exitOnFailure, runWithStdout
from checkm.hmmerAligner import HmmerAligner
from checkm.binStatistics import BinStatistics


def haveProgram(program):
    return shutil.which(program) is not None


class Hit():
    def __init__(self, target_name, full_e_value, full_score):
        self.target_name = target_name
        self.full_e_value = full_e_value
        self.full_score = full_score


class VerifyReturnCodes(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger('timestamp')
        self.tmpDir = tempfile.mkdtemp(prefix='checkm test ')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testReturnCodeStr(self):
        """Verify description of subprocess return codes."""
        self.assertEqual(returnCodeStr(1), 'exit code 1')
        self.assertEqual(returnCodeStr(-9), 'signal SIGKILL')

    def testExitOnFailure(self):
        """Verify failures exit with status 1, regardless of the program's exit code or signal."""
        exitOnFailure(self.logger, 'prog', 0)

        for returnCode in [1, 2, 255, 256, -9]:
            with self.assertLogs('timestamp', level='ERROR'):
                with self.assertRaises(SystemExit) as cm:
                    exitOnFailure(self.logger, 'prog', returnCode)
            self.assertEqual(cm.exception.code, 1)

    def testRunWithStdout(self):
        """Verify exit codes are reported directly (not as wait statuses) and output is written to paths with spaces."""
        outputFile = os.path.join(self.tmpDir, 'output file.txt')

        self.assertEqual(runWithStdout(['sh', '-c', 'echo first'], outputFile), 0)
        self.assertEqual(runWithStdout(['sh', '-c', 'echo second'], outputFile, append=True), 0)
        with open(outputFile) as f:
            self.assertEqual(f.read(), 'first\nsecond\n')

        self.assertEqual(runWithStdout(['sh', '-c', 'exit 1'], outputFile), 1)
        self.assertEqual(runWithStdout(['sh', '-c', 'kill -9 $$'], outputFile), -9)

    @unittest.skipUnless(haveProgram('hmmfetch') and haveProgram('hmmsearch'), 'HMMER is not installed')
    def testHmmerFailureExitCode(self):
        """Verify a failing HMMER program exits with status 1."""
        HF = HMMERRunner(mode='fetch')
        with self.assertLogs('timestamp', level='ERROR'):
            with self.assertRaises(SystemExit) as cm:
                HF.fetch(os.path.join(self.tmpDir, 'missing models.hmm'), 'PF00001', os.path.join(self.tmpDir, 'out.hmm'))
        self.assertEqual(cm.exception.code, 1)

    @unittest.skipUnless(haveProgram('hmmfetch') and haveProgram('hmmsearch'), 'HMMER is not installed')
    def testHmmerPathsWithSpaces(self):
        """Verify HMMER programs can be run on paths containing spaces."""
        if not os.path.exists(DefaultValues.PHYLO_HMM_MODELS):
            self.skipTest('CheckM reference data is not available')

        hmmFile = os.path.join(self.tmpDir, 'models with spaces.hmm')
        shutil.copyfile(DefaultValues.PHYLO_HMM_MODELS, hmmFile)

        HF = HMMERRunner(mode='fetch')
        HF.index(hmmFile)
        self.assertTrue(os.path.exists(hmmFile + '.ssi'))

        modelFile = os.path.join(self.tmpDir, 'one model.hmm')
        HF.fetch(hmmFile, 'PF00164.20', modelFile)
        with open(modelFile) as f:
            self.assertIn('PF00164.20', f.read())


    @unittest.skipUnless(haveProgram('hmmfetch') and haveProgram('hmmsearch'), 'HMMER is not installed')
    def testExtractModelsFailure(self):
        """Verify a failure extracting HMMs for alignment exits with status 1."""
        if not os.path.exists(DefaultValues.PHYLO_HMM_MODELS):
            self.skipTest('CheckM reference data is not available')

        HA = HmmerAligner(threads=2)
        hmmModelFiles = {}
        with self.assertLogs('timestamp', level='ERROR'):
            with self.assertRaises(SystemExit) as cm:
                HA._makeAlignmentModels(DefaultValues.PHYLO_HMM_MODELS, ['PF00164.20', 'NOT_A_MODEL'], hmmModelFiles, bReportProgress=False)
        self.assertEqual(cm.exception.code, 1)

        HA._removeFiles(hmmModelFiles.values())

    @unittest.skipUnless(haveProgram('hmmfetch') and haveProgram('hmmsearch'), 'HMMER is not installed')
    def testWorkerFailureDoesNotHang(self):
        """Verify the run exits, rather than hanging, when workers fail with many unprocessed items queued."""
        script = (
            "import logging, os, tempfile\n"
            "from checkm.hmmerAligner import HmmerAligner\n"
            "logging.getLogger('timestamp').addHandler(logging.NullHandler())\n"
            "hmmModelFiles = {}\n"
            "try:\n"
            "    HmmerAligner(threads=2)._makeAlignmentModels(os.path.join(tempfile.gettempdir(), 'missing.hmm'),\n"
            "                                                 ['MODEL_%d' % i for i in range(5000)], hmmModelFiles, False)\n"
            "finally:\n"
            "    HmmerAligner(threads=1)._removeFiles(hmmModelFiles.values())\n"
        )

        env = dict(os.environ)
        env['PYTHONPATH'] = os.pathsep.join([os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))] +
                                            ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))
        env['TMPDIR'] = self.tmpDir
        proc = subprocess.run([sys.executable, '-c', script], env=env, timeout=300,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.assertEqual(proc.returncode, 1)


class VerifyTopHit(unittest.TestCase):
    def testTopHit(self):
        """Verify the best hit (lowest e-value, then highest score) is selected for --ali."""
        HA = HmmerAligner(threads=1)

        hits = [Hit('worst', 1e-5, 20.0), Hit('best', 1e-50, 150.0), Hit('middle', 1e-20, 80.0)]
        self.assertEqual(HA._topHit(hits).target_name, 'best')

        # ties in e-value are broken by the highest bit score
        hits = [Hit('a', 0.0, 900.0), Hit('b', 0.0, 1000.0), Hit('c', 1e-10, 50.0)]
        self.assertEqual(HA._topHit(hits).target_name, 'b')

        # hit list is not reordered
        self.assertEqual([h.target_name for h in hits], ['a', 'b', 'c'])


class VerifyEmptyBinStatistics(unittest.TestCase):
    def testEmptyBins(self):
        """Verify sequence statistics for bins with no sequences or only ambiguous bases."""
        binStats = BinStatistics(threads=1)

        stats = binStats.calculateSeqStats({})
        self.assertEqual(stats, (0, 0, 0, 0, 0, 0.0, 0.0, 0, 0))

        seqStats = defaultdict(dict)
        stats = binStats.calculateSeqStats({'S1': 'NNNNNNNNNN', 'S2': 'NNNNN'}, seqStats)
        (maxScaffoldLen, maxContigLen, genomeSize, scaffoldN50, contigN50,
         scaffoldAvgLen, contigAvgLen, numContigs, numAmbiguousBases) = stats
        self.assertEqual(maxScaffoldLen, 10)
        self.assertEqual(maxContigLen, 0)
        self.assertEqual(genomeSize, 15)
        self.assertEqual(scaffoldN50, 10)
        self.assertEqual(contigN50, 0)
        self.assertAlmostEqual(scaffoldAvgLen, 7.5)
        self.assertEqual(contigAvgLen, 0.0)
        self.assertEqual(numContigs, 0)
        self.assertEqual(numAmbiguousBases, 15)
        self.assertEqual(seqStats['S1']['# contigs'], 0)

        self.assertEqual(binStats.calculateGC({}), (0.0, 0.0))

    def testCodingDensityEmptyBin(self):
        """Verify coding density of a bin with no bases is 0 rather than a division by zero."""
        binStats = BinStatistics(threads=1)

        tmpDir = tempfile.mkdtemp()
        try:
            open(os.path.join(tmpDir, DefaultValues.PRODIGAL_GFF), 'w').close()
            open(os.path.join(tmpDir, DefaultValues.PRODIGAL_AA), 'w').close()

            codingDensity, _translationTable, numGenes = binStats.calculateCodingDensity(tmpDir, {}, 0)
            self.assertEqual(codingDensity, 0.0)
            self.assertEqual(numGenes, 0)
        finally:
            shutil.rmtree(tmpDir)


@unittest.skipUnless(haveProgram('prodigal') and haveProgram('hmmsearch'), 'prodigal or HMMER is not installed')
class VerifyGenelessBins(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm test ')

        # 200 bp sequence in which prodigal finds no genes
        self.tinyBin = os.path.join(self.tmpDir, 'tiny bin.fna')
        with open(self.tinyBin, 'w') as f:
            f.write('>tiny\n' + 'ACGTTGCA' * 25 + '\n')

        self.emptyBin = os.path.join(self.tmpDir, 'empty bin.fna')
        open(self.emptyBin, 'w').close()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testProdigalNoGenes(self):
        """Verify prodigal finding no genes, or being given no sequences, is not treated as a failure."""
        from checkm.prodigal import ProdigalRunner

        for binFile, outDir in [(self.tinyBin, 'tiny out'), (self.emptyBin, 'empty out')]:
            outDir = os.path.join(self.tmpDir, outDir)
            os.makedirs(outDir)

            prodigal = ProdigalRunner(outDir)
            prodigal.run(binFile, bNucORFs=True)

            self.assertTrue(os.path.exists(prodigal.aaGeneFile))
            self.assertEqual(os.path.getsize(prodigal.aaGeneFile), 0)
            self.assertTrue(os.path.exists(prodigal.gffFile))
            self.assertFalse(prodigal.areORFsCalled(bNucORFs=False))

    def testProdigalEmptyRecord(self):
        """Verify a record without sequence data doesn't make prodigal fail."""
        from checkm.prodigal import ProdigalRunner

        # prodigal exits with 'no input sequences to analyze' if the first record is empty
        binFile = os.path.join(self.tmpDir, 'empty record.fna')
        with open(binFile, 'w') as f:
            f.write('>empty\n>seq\n' + 'ACGTTGCA' * 25 + '\n')

        outDir = os.path.join(self.tmpDir, 'empty record out')
        os.makedirs(outDir)

        prodigal = ProdigalRunner(outDir)
        with self.assertLogs('timestamp', level='WARNING') as cm:
            prodigal.run(binFile, bNucORFs=True)

        self.assertTrue(any('without sequence data' in msg for msg in cm.output))
        self.assertTrue(os.path.exists(prodigal.aaGeneFile))
        self.assertTrue(os.path.exists(prodigal.gffFile))

    def testMarkerGeneFinderSkipsGenelessBins(self):
        """Verify bins without genes are skipped with a warning and other bins are processed."""
        from checkm.markerGeneFinder import MarkerGeneFinder

        if not os.path.exists(DefaultValues.PHYLO_HMM_MODELS):
            self.skipTest('CheckM reference data is not available')

        # called genes: a bin with a few genes, and bins without any genes
        geneBins = []
        for binId, seqs in [('b_bin', {'gene_1': 'MKVLAAGIVGLLLAPTAQAQ' * 5}),
                            ('a_bin', {}),
                            ('c_bin', {'gene_1': 'MSTNPKPQRKTKRNTNRRPQ' * 5, 'gene_2': 'MADEEKLPPGWEKRMSRSSG' * 5})]:
            binFile = os.path.join(self.tmpDir, binId + '.faa')
            with open(binFile, 'w') as f:
                for seqId, seq in seqs.items():
                    f.write('>%s\n%s\n' % (seqId, seq))
            geneBins.append(binFile)

        outDir = os.path.join(self.tmpDir, 'out dir')
        os.makedirs(os.path.join(outDir, 'bins'))

        mgf = MarkerGeneFinder(threads=2)
        with self.assertLogs('timestamp', level='WARNING') as cm:
            binIdToModels = mgf.find(geneBins, outDir,
                                     DefaultValues.HMMER_TABLE_PHYLO_OUT,
                                     DefaultValues.HMMER_PHYLO_OUT,
                                     DefaultValues.PHYLO_HMM_MODELS,
                                     False, False, True)

        self.assertEqual(list(binIdToModels.keys()), ['b_bin', 'c_bin'])
        self.assertEqual(len([msg for msg in cm.output if 'Skipping bin a_bin' in msg]), 1)
        for binId in ['b_bin', 'c_bin']:
            self.assertTrue(os.path.exists(os.path.join(outDir, 'bins', binId, DefaultValues.HMMER_TABLE_PHYLO_OUT)))

        # gene calling: bins in which prodigal finds no genes are skipped
        outDir = os.path.join(self.tmpDir, 'out dir nt')
        os.makedirs(os.path.join(outDir, 'bins'))
        with self.assertLogs('timestamp', level='WARNING') as cm:
            binIdToModels = mgf.find([self.tinyBin, self.emptyBin], outDir,
                                     DefaultValues.HMMER_TABLE_PHYLO_OUT,
                                     DefaultValues.HMMER_PHYLO_OUT,
                                     DefaultValues.PHYLO_HMM_MODELS,
                                     False, False, False)

        self.assertEqual(binIdToModels, {})
        self.assertEqual(len([msg for msg in cm.output if 'no genes were identified by Prodigal' in msg]), 2)

    def testMarkerGeneFinderWorkerFailure(self):
        """Verify a failed worker aborts the run with exit status 1 rather than silently dropping the bin."""
        from checkm.markerGeneFinder import MarkerGeneFinder

        # protein FASTA that hmmsearch cannot parse
        badBin = os.path.join(self.tmpDir, 'bad.faa')
        with open(badBin, 'w') as f:
            f.write('this is not a FASTA file\n')

        outDir = os.path.join(self.tmpDir, 'out dir')
        os.makedirs(os.path.join(outDir, 'bins'))

        if not os.path.exists(DefaultValues.PHYLO_HMM_MODELS):
            self.skipTest('CheckM reference data is not available')

        mgf = MarkerGeneFinder(threads=1)
        with self.assertLogs('timestamp', level='ERROR'):
            with self.assertRaises(SystemExit) as cm:
                mgf.find([badBin], outDir,
                         DefaultValues.HMMER_TABLE_PHYLO_OUT,
                         DefaultValues.HMMER_PHYLO_OUT,
                         DefaultValues.PHYLO_HMM_MODELS,
                         False, False, True)
        self.assertEqual(cm.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
