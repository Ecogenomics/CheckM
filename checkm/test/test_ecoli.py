###############################################################################
#
# test_ecoli.py - process E.coli K12-W3310 genome to verify operation of CheckM
#
# Note: This test must be initiated using checkm test <output dir>.
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
import logging
import shutil

from checkm.resultsParser import ResultsParser
from checkm.markerSets import MarkerSet
from checkm.defaultValues import DefaultValues
from checkm.common import makeSurePathExists, checkFileExists

import numpy as np


class Options():
    pass


class VerifyEcoli():
    def __init__(self):
        self.logger = logging.getLogger('timestamp')

    def run(self, parser, outputDir):
        """Run standard E. coli genome to verify operation of CheckM."""

        ecoliFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'test_data', '637000110.fna')
        checkFileExists(ecoliFile)

        options = Options()
        options.threads = 1
        options.pplacer_threads = 1
        options.extension = 'fna'
        options.bQuiet = True
        options.output_dir = os.path.join(outputDir, 'results')
        if os.path.exists(options.output_dir):
            shutil.rmtree(options.output_dir)
        makeSurePathExists(options.output_dir)

        self.logger.info('[Step 1]: Verifying tree command.')
        options.bKeepAlignment = False
        options.bNucORFs = False
        options.bCalledGenes = False
        options.bReducedTree = True
        options.bin_input = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'test_data')
        parser.tree(options)
        self.verifyTree(options.output_dir)
        self.logger.info('[Passed]')

        self.logger.info('[Step 2]: Verifying tree_qa command.')
        options.tree_dir = options.output_dir
        options.out_format = 1
        options.file = os.path.join(options.output_dir, 'tree_qa_test.tsv')
        options.bTabTable = True
        parser.treeQA(options)
        self.verifyTreeQA(options.file)
        self.logger.info('[Passed]')

        self.logger.info('[Step 3]: Verifying lineage_set command.')
        options.marker_file = os.path.join(options.output_dir, 'lineage_set_test.tsv')
        options.bForceDomain = False
        options.bootstrap = 0
        options.num_genomes_markers = 30
        options.num_genomes_refine = 5
        options.bNoLineageSpecificRefinement = False

        options.bRequireTaxonomy = False
        options.unique = 10
        options.multi = 10
        parser.lineageSet(options)
        self.verifyLineageSet(options.marker_file, options.bRequireTaxonomy)

        options.bRequireTaxonomy = True
        parser.lineageSet(options)
        self.verifyLineageSet(options.marker_file, options.bRequireTaxonomy)
        self.logger.info('[Passed]')

        self.logger.info('[Step 4]: Verifying analyze command.')
        options.bAlignTopHit = False
        parser.analyze(options)
        self.verifyAnalyze(options.output_dir)
        self.logger.info('[Passed]')

        self.logger.info('[Step 5]: Verifying qa command.')
        options.alignment_file = None
        options.analyze_dir = options.output_dir
        options.out_format = 1
        options.exclude_markers = None
        options.bSkipPseudoGeneCorrection = False
        options.bSkipAdjCorrection = False
        options.file = os.path.join(options.output_dir, 'qa_test.tsv')
        options.bIndividualMarkers = False
        options.bIgnoreThresholds = False
        options.aai_strain = 0.9
        options.e_value = 1e-10
        options.length = 0.7
        options.coverage_file = None
        options.bTabTable = True
        parser.qa(options)
        self.verifyQA(options.file)
        self.logger.info('[Passed]')

    def verifyTree(self, outdir):
        """Verify output of tree command."""

        # verify bin stats using independently verified ground truth values
        binStats = ResultsParser(None).parseBinStats(outdir, DefaultValues.BIN_STATS_PHYLO_OUT)

        np.testing.assert_almost_equal(binStats['637000110']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(binStats['637000110']['GC std'], 0.0, err_msg="Failed GC std test")
        # np.testing.assert_almost_equal(binStats['637000110']['Coding density'], 0.8775, decimal=3, err_msg="Failed coding density test") # depends on exact version of prodigal
        np.testing.assert_almost_equal(binStats['637000110']['# contigs'], 1, err_msg="Failed # contigs test")
        np.testing.assert_almost_equal(binStats['637000110']['# scaffolds'], 1, err_msg="Failed # scaffolds test")
        np.testing.assert_equal(binStats['637000110']['Longest contig'], 4646332, err_msg="Failed longest contig test")
        np.testing.assert_equal(binStats['637000110']['Longest scaffold'], 4646332, err_msg="Failed longest scaffold test")
        # np.testing.assert_equal(binStats['637000110']['# predicted genes'], 4327, err_msg="Failed # predicted genes test") # depends on exact version of prodigal
        np.testing.assert_equal(binStats['637000110']['N50 (contigs)'], 4646332, err_msg="Failed N50 (contigs) test")
        np.testing.assert_equal(binStats['637000110']['N50 (scaffolds)'], 4646332, err_msg="Failed N50 (scaffolds) test")
        np.testing.assert_equal(binStats['637000110']['Genome size'], 4646332, err_msg="Failed genome size test")

        # verify sequence stats using  independently verified ground truth values
        # [The sequence stats file is not generated in CheckM >= v0.9.8 in order to reduce memory requirements.]
        # with open(os.path.join(outdir, 'storage', DefaultValues.SEQ_STATS_PHYLO_OUT), 'r') as f:
        #    s = f.read()
        #    seqStats = ast.literal_eval(s)

        # np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['Total contig length'], 4646332, err_msg="Failed total contig length test")
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['Coding bases'], 4077069, err_msg="Failed coding bases test") # depends on exact version of prodigal
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['# ORFs'], 4326, err_msg="Failed # genes test") # depends on exact version of prodigal
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['Length'], 4646332, err_msg="Failed length test")
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['# contigs'], 1, err_msg="Failed # contigs test")

    def verifyTreeQA(self, qaTableFile):
        """Verify output of tree QA command."""

        with open(qaTableFile) as f:
            f.readline()  # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')

        np.testing.assert_almost_equal(int(lineSplit[0]), 637000110, err_msg="Failed genome ID test")
        np.testing.assert_almost_equal(int(lineSplit[1]), 43, err_msg="Failed # markers test")

        family = None
        if len(lineSplit) >= 4:
            taxonomy = lineSplit[3].split(';')
            if len(taxonomy) >= 5:
                family = lineSplit[3].split(';')[4].strip()
        assert(family == 'f__Enterobacteriaceae')

    def verifyLineageSet(self, markerSetFile, bRequireTaxonomy):
        """Verify output of lineage set command."""

        with open(markerSetFile) as f:
            f.readline()  # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')
                    binId = lineSplit[0]
                    _numMarkers = int(lineSplit[1])
                    uid = lineSplit[2]
                    lineage = lineSplit[3]
                    numGenomes = int(lineSplit[4])
                    _markerSet = MarkerSet(uid, lineage, numGenomes, eval(lineSplit[5].rstrip()))

        np.testing.assert_almost_equal(int(binId), 637000110, err_msg="Failed bin ID test")
        if not bRequireTaxonomy:
            # this is unstable as it depends on HMMER and prodigal
            # np.testing.assert_equal(markerSet.numSets(), 336, err_msg="Failed # marker set test")
            # np.testing.assert_equal(markerSet.numMarkers(), 1173, err_msg="Failed # markers test")
            pass
        else:
            # np.testing.assert_equal(markerSet.numSets(), 282, err_msg="Failed # marker set test")
            # np.testing.assert_equal(markerSet.numMarkers(), 1254, err_msg="Failed # markers test")
            pass

    def verifyAnalyze(self, outdir):
        """Verify output of analyze command."""

        # verify bin stats using independently verified ground truth values
        binStats = ResultsParser(None).parseBinStats(outdir, DefaultValues.BIN_STATS_OUT)

        np.testing.assert_almost_equal(binStats['637000110']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(binStats['637000110']['GC std'], 0.0, err_msg="Failed GC std test")
        # np.testing.assert_almost_equal(binStats['637000110']['Coding density'], 0.877, decimal=3, err_msg="Failed coding density test") # depends on exact version of prodigal
        np.testing.assert_equal(binStats['637000110']['# contigs'], 1, err_msg="Failed # contigs test")
        np.testing.assert_equal(binStats['637000110']['# scaffolds'], 1, err_msg="Failed # scaffolds test")
        np.testing.assert_equal(binStats['637000110']['Longest contig'], 4646332, err_msg="Failed longest contig test")
        np.testing.assert_equal(binStats['637000110']['Longest scaffold'], 4646332, err_msg="Failed longest scaffold test")
        # np.testing.assert_equal(binStats['637000110']['# predicted genes'], 4326, err_msg="Failed # predicted genes test") # depends on exact version of prodigal
        np.testing.assert_equal(binStats['637000110']['N50 (contigs)'], 4646332, err_msg="Failed N50 (contigs) test")
        np.testing.assert_equal(binStats['637000110']['N50 (scaffolds)'], 4646332, err_msg="Failed N50 (scaffolds) test")
        np.testing.assert_equal(binStats['637000110']['Genome size'], 4646332, err_msg="Failed genome size test")

        # verify sequence stats using  independently verified ground truth values
        # [The sequence stats file is not generated in CheckM >= v0.9.8 in order to reduce memory requirements.]
        # with open(os.path.join(outdir, 'storage', DefaultValues.SEQ_STATS_OUT), 'r') as f:
        #    s = f.read()
        #    seqStats = ast.literal_eval(s)

        # np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['Total contig length'], 4646332, err_msg="Failed total contig length test")
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['Coding bases'], 4077069, err_msg="Failed coding bases test") # depends on exact version of prodigal
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['# ORFs'], 4326, err_msg="Failed # genes test") # depends on exact version of prodigal
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['Length'], 4646332, err_msg="Failed length test")
        # np.testing.assert_equal(seqStats['637000110']['AC_000091']['# contigs'], 1, err_msg="Failed # contigs test")

    def verifyQA(self, qaTableFile):
        """Verify output of qa command."""

        with open(qaTableFile) as f:
            f.readline()  # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')

        np.testing.assert_equal(int(lineSplit[0]), 637000110, err_msg="Failed genome ID test")
        # np.testing.assert_equal(lineSplit[1], 'f__Enterobacteriaceae', err_msg="Failed lineage test")    # depend on exact version of prodigal
        # np.testing.assert_equal(int(lineSplit[2]), 134, err_msg="Failed # genomes")
        # np.testing.assert_equal(int(lineSplit[3]), 1173, err_msg="Failed # markers test")
        # np.testing.assert_almost_equal(int(lineSplit[4]), 336, err_msg="Failed # marker sets test")
        # np.testing.assert_almost_equal(float(lineSplit[11]), 99.98, decimal=2, err_msg="Failed completeness test")
        # np.testing.assert_almost_equal(float(lineSplit[12]), 0.04, decimal=2, err_msg="Failed contamination test")