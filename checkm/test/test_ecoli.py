###############################################################################
#
# test_ecoli.py - process E.coli K12-W3310 genome to verify operation of CheckM
#
# Note: This test must be initiated separately from 'nosetests'.
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
import ast
import shutil

from checkm.markerSets import MarkerSet
from checkm.defaultValues import DefaultValues
from checkm.common import makeSurePathExists, checkFileExists
from checkm.main import OptionsParser

import numpy as np

class VerifyEcoli():
    def __init__(self):
        pass

    def verifyTree(self, outdir):
        """Verify output of tree command."""

        # verify bin stats using independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.BIN_STATS_PHYLO_OUT), 'r') as f:
            s = f.read()
            binStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(binStats['637000110']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(binStats['637000110']['GC std'], 0.0, err_msg="Failed GC std test")
        np.testing.assert_almost_equal(binStats['637000110']['Coding density'], 0.877, decimal=3, err_msg="Failed coding density test")
        np.testing.assert_almost_equal(binStats['637000110']['# contigs'], 1, err_msg="Failed # contigs test")
        np.testing.assert_almost_equal(binStats['637000110']['# scaffolds'], 1, err_msg="Failed # scaffolds test")
        np.testing.assert_equal(binStats['637000110']['Longest contig'], 4646332, err_msg="Failed longest contig test")
        np.testing.assert_equal(binStats['637000110']['Longest scaffold'], 4646332, err_msg="Failed longest scaffold test")
        np.testing.assert_equal(binStats['637000110']['# predicted ORFs'], 4326, err_msg="Failed # predicted ORFs test")
        np.testing.assert_equal(binStats['637000110']['N50 (contigs)'], 4646332, err_msg="Failed N50 (contigs) test")
        np.testing.assert_equal(binStats['637000110']['N50 (scaffolds)'], 4646332, err_msg="Failed N50 (scaffolds) test")
        np.testing.assert_equal(binStats['637000110']['Genome size'], 4646332, err_msg="Failed genome size test")

        # verify sequence stats using  independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.SEQ_STATS_PHYLO_OUT), 'r') as f:
            s = f.read()
            seqStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['Total contig length'], 4646332, err_msg="Failed total contig length test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['Coding bases'], 4077069, err_msg="Failed coding bases test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['# ORFs'], 4326, err_msg="Failed # ORFs test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['Length'], 4646332, err_msg="Failed length test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['# contigs'], 1, err_msg="Failed # contigs test")

    def verifyTreeQA(self, qaTableFile):
        """Verify output of tree QA command."""

        with open(qaTableFile) as f:
            f.readline() # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')

        np.testing.assert_almost_equal(int(lineSplit[0]), 637000110, err_msg="Failed genome ID test")
        np.testing.assert_almost_equal(int(lineSplit[1]), 43, err_msg="Failed # markers test")

        family = lineSplit[2].split(';')[4].strip()
        assert(family == 'f__Enterobacteriaceae')

    def verifyLineageSet(self, markerSetFile, bRequireTaxonomy):
        """Verify output of lineage set command."""

        with open(markerSetFile) as f:
            f.readline() # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')
                    binId = lineSplit[0]
                    numMarkers = int(lineSplit[1])
                    uid = lineSplit[2]
                    lineage = lineSplit[3]
                    numGenomes = int(lineSplit[4])
                    markerSet = MarkerSet(uid, lineage, numGenomes, eval(lineSplit[5].rstrip()))

        np.testing.assert_almost_equal(int(binId), 637000110, err_msg="Failed bin ID test")
        if not bRequireTaxonomy:
            np.testing.assert_equal(markerSet.numSets(), 245, err_msg="Failed # marker set test")
            np.testing.assert_equal(markerSet.numMarkers(), 1945, err_msg="Failed # markers test")
            assert(uid == 'UID5199')
        else:
            np.testing.assert_equal(markerSet.numSets(), 282, err_msg="Failed # marker set test")
            np.testing.assert_equal(markerSet.numMarkers(), 1254, err_msg="Failed # markers test")
            assert(lineage == 'f__Enterobacteriaceae')

    def verifyAnalyze(self, outdir):
        """Verify output of analyze command."""

                # verify bin stats using independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.BIN_STATS_OUT), 'r') as f:
            s = f.read()
            binStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(binStats['637000110']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(binStats['637000110']['GC std'], 0.0, err_msg="Failed GC std test")
        np.testing.assert_almost_equal(binStats['637000110']['Coding density'], 0.877, decimal=3, err_msg="Failed coding density test")
        np.testing.assert_equal(binStats['637000110']['# contigs'], 1, err_msg="Failed # contigs test")
        np.testing.assert_equal(binStats['637000110']['# scaffolds'], 1, err_msg="Failed # scaffolds test")
        np.testing.assert_equal(binStats['637000110']['Longest contig'], 4646332, err_msg="Failed longest contig test")
        np.testing.assert_equal(binStats['637000110']['Longest scaffold'], 4646332, err_msg="Failed longest scaffold test")
        np.testing.assert_equal(binStats['637000110']['# predicted ORFs'], 4326, err_msg="Failed # predicted ORFs test")
        np.testing.assert_equal(binStats['637000110']['N50 (contigs)'], 4646332, err_msg="Failed N50 (contigs) test")
        np.testing.assert_equal(binStats['637000110']['N50 (scaffolds)'], 4646332, err_msg="Failed N50 (scaffolds) test")
        np.testing.assert_equal(binStats['637000110']['Genome size'], 4646332, err_msg="Failed genome size test")

        # verify sequence stats using  independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.SEQ_STATS_OUT), 'r') as f:
            s = f.read()
            seqStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['Total contig length'], 4646332, err_msg="Failed total contig length test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['Coding bases'], 4077069, err_msg="Failed coding bases test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['# ORFs'], 4326, err_msg="Failed # ORFs test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['Length'], 4646332, err_msg="Failed length test")
        np.testing.assert_equal(seqStats['637000110']['AC_000091']['# contigs'], 1, err_msg="Failed # contigs test")

    def verifyQA(self, qaTableFile):
        """Verify output of qa command."""

        with open(qaTableFile) as f:
            f.readline() # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')

        np.testing.assert_equal(int(lineSplit[0]), 637000110, err_msg="Failed genome ID test")
        np.testing.assert_equal(lineSplit[1], 'f__Enterobacteriaceae', err_msg="Failed lineage test")
        np.testing.assert_equal(int(lineSplit[2]), 134, err_msg="Failed # genomes")
        np.testing.assert_equal(int(lineSplit[3]), 1100, err_msg="Failed # markers test")
        np.testing.assert_almost_equal(int(lineSplit[4]), 327, err_msg="Failed # marker sets test")
        np.testing.assert_almost_equal(float(lineSplit[11]), 99.98, decimal=2, err_msg="Failed completeness test")
        np.testing.assert_almost_equal(float(lineSplit[12]), 0.04, decimal=2, err_msg="Failed contamination test")
        
class Options():
    pass
        
if __name__ == "__main__":
    """Bin compare command"""
    print ''
    print '*******************************************************************************'
    print '[CheckM - test] Processing E.coli K12-W3310 and verify operation of CheckM.'
    print '*******************************************************************************'

    ecoliFile = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'test_data', '637000110.fna')
    checkFileExists(ecoliFile)

    parser = OptionsParser()
    verify = VerifyEcoli()
    
    options = Options()
    options.threads = 1
    options.extension = 'fna'
    options.bQuiet = True
    options.out_folder = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'test_data', 'results')
    if os.path.exists(options.out_folder):
        shutil.rmtree(options.out_folder)
    makeSurePathExists(options.out_folder)

    print '  [Step 1]: Verifying tree command.'
    options.bKeepAlignment = False
    options.bNucORFs = False
    options.bin_folder = os.path.join(DefaultValues.CHECKM_DATA_DIR, 'test_data')
    parser.tree(options)
    verify.verifyTree(options.out_folder)
    print '    Passed.'

    print '  [Step 2]: Verifying tree_qa command.'
    options.tree_folder = options.out_folder
    options.out_format = 1
    options.file = os.path.join(options.out_folder, 'tree_qa_test.tsv')
    options.bTabTable = True
    parser.treeQA(options)
    verify.verifyTreeQA(options.file)
    print '    Passed.'

    print '  [Step 3]: Verifying lineage_set command.'
    options.marker_file = os.path.join(options.out_folder, 'lineage_set_test.tsv')
    options.bForceDomain = False
    options.bootstrap = 0
    options.num_genomes_markers = 30
    options.num_genomes_refine = 5
    options.bNoLineageSpecificRefinement = False

    options.bRequireTaxonomy = False
    parser.lineageSet(options)
    verify.verifyLineageSet(options.marker_file, options.bRequireTaxonomy)

    options.bRequireTaxonomy = True
    parser.lineageSet(options)
    verify.verifyLineageSet(options.marker_file, options.bRequireTaxonomy)
    print '    Passed.'

    print '  [Step 4]: Verifying analyze command.'
    options.bAlignTopHit = False    
    parser.analyze(options)
    verify.verifyAnalyze(options.out_folder)
    print '    Passed.'

    print '  [Step 5]: Verifying qa command.'
    options.alignment_file = None
    options.analyze_folder = options.out_folder
    options.out_format = 1
    options.file = os.path.join(options.out_folder, 'qa_test.tsv')
    options.bIndividualMarkers = False
    options.bSkipOrfCorrection = False
    options.bIgnoreThresholds = False
    options.aai_strain = 0.95
    options.e_value = 1e-10
    options.length = 0.7
    options.coverage_file = None
    options.bTabTable = True
    parser.qa(options)
    verify.verifyQA(options.file)
    print '    Passed.'
  
