###############################################################################
#
# testing.py - test operation of CheckM
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

from checkm.markerSets import MarkerSet
from checkm.defaultValues import DefaultValues

import numpy as np

class Testing():
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
        np.testing.assert_almost_equal(binStats['637000110']['Longest contig'], 4646332, err_msg="Failed longest contig test")
        np.testing.assert_almost_equal(binStats['637000110']['Longest scaffold'], 4646332, err_msg="Failed longest scaffold test")
        np.testing.assert_almost_equal(binStats['637000110']['# predicted ORFs'], 4326, err_msg="Failed # predicted ORFs test")
        np.testing.assert_almost_equal(binStats['637000110']['N50 (contigs)'], 4646332, err_msg="Failed N50 (contigs) test")
        np.testing.assert_almost_equal(binStats['637000110']['N50 (scaffolds)'], 4646332, err_msg="Failed N50 (scaffolds) test")
        np.testing.assert_almost_equal(binStats['637000110']['Genome size'], 4646332, err_msg="Failed genome size test")

        # verify sequence stats using  independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.SEQ_STATS_PHYLO_OUT), 'r') as f:
            s = f.read()
            seqStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['Total contig length'], 4646332, err_msg="Failed total contig length test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['Coding bases'], 4077069, err_msg="Failed coding bases test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['# ORFs'], 4326, err_msg="Failed # ORFs test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['Length'], 4646332, err_msg="Failed length test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['# contigs'], 1, err_msg="Failed # contigs test")

    def verifyTreeQA(self, qaTableFile):
        """Verify output of tree QA command."""

        with open(qaTableFile) as f:
            f.readline() # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')

        np.testing.assert_almost_equal(int(lineSplit[0]), 637000110, err_msg="Failed genome ID test")
        np.testing.assert_almost_equal(int(lineSplit[1]), 36, err_msg="Failed # markers test")

        family = lineSplit[2].split(';')[4]
        assert(family == 'f__Enterobacteriaceae')

    def verifyLineageSet(self, markerSetFile, bRequireTaxonomy):
        """Verify output of lineage set command."""

        with open(markerSetFile) as f:
            f.readline() # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')
                    binId = lineSplit[0]
                    binPlacement = lineSplit[1]
                    markerSet = MarkerSet(eval(lineSplit[4].rstrip()))

        np.testing.assert_almost_equal(int(binId), 637000110, err_msg="Failed bin ID test")
        if not bRequireTaxonomy:
            np.testing.assert_almost_equal(markerSet.numSets(), 286, err_msg="Failed # marker set test")
            np.testing.assert_almost_equal(markerSet.numMarkers(), 1290, err_msg="Failed # markers test")
            assert(binPlacement == 'UID3189')
        else:
            np.testing.assert_almost_equal(markerSet.numSets(), 282, err_msg="Failed # marker set test")
            np.testing.assert_almost_equal(markerSet.numMarkers(), 1254, err_msg="Failed # markers test")
            assert(binPlacement == 'g__Shigella')

    def verifyAnalyze(self, outdir):
        """Verify output of analyze command."""

                # verify bin stats using independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.BIN_STATS_OUT), 'r') as f:
            s = f.read()
            binStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(binStats['637000110']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(binStats['637000110']['GC std'], 0.0, err_msg="Failed GC std test")
        np.testing.assert_almost_equal(binStats['637000110']['Coding density'], 0.877, decimal=3, err_msg="Failed coding density test")
        np.testing.assert_almost_equal(binStats['637000110']['# contigs'], 1, err_msg="Failed # contigs test")
        np.testing.assert_almost_equal(binStats['637000110']['# scaffolds'], 1, err_msg="Failed # scaffolds test")
        np.testing.assert_almost_equal(binStats['637000110']['Longest contig'], 4646332, err_msg="Failed longest contig test")
        np.testing.assert_almost_equal(binStats['637000110']['Longest scaffold'], 4646332, err_msg="Failed longest scaffold test")
        np.testing.assert_almost_equal(binStats['637000110']['# predicted ORFs'], 4326, err_msg="Failed # predicted ORFs test")
        np.testing.assert_almost_equal(binStats['637000110']['N50 (contigs)'], 4646332, err_msg="Failed N50 (contigs) test")
        np.testing.assert_almost_equal(binStats['637000110']['N50 (scaffolds)'], 4646332, err_msg="Failed N50 (scaffolds) test")
        np.testing.assert_almost_equal(binStats['637000110']['Genome size'], 4646332, err_msg="Failed genome size test")

        # verify sequence stats using  independently verified ground truth values
        with open(os.path.join(outdir, 'storage', DefaultValues.SEQ_STATS_OUT), 'r') as f:
            s = f.read()
            seqStats = ast.literal_eval(s)

        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['GC'], 0.508, decimal=3, err_msg="Failed GC test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['Total contig length'], 4646332, err_msg="Failed total contig length test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['Coding bases'], 4077069, err_msg="Failed coding bases test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['# ORFs'], 4326, err_msg="Failed # ORFs test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['Length'], 4646332, err_msg="Failed length test")
        np.testing.assert_almost_equal(seqStats['637000110']['AC_000091']['# contigs'], 1, err_msg="Failed # contigs test")

    def verifyQA(self, qaTableFile):
        """Verify output of qa command."""

        with open(qaTableFile) as f:
            f.readline() # skip header

            for line in f:
                if line.strip() != '':
                    lineSplit = line.split('\t')

        np.testing.assert_almost_equal(int(lineSplit[0]), 637000110, err_msg="Failed genome ID test")
        np.testing.assert_almost_equal(int(lineSplit[1]), 1254, err_msg="Failed # markers test")
        np.testing.assert_almost_equal(int(lineSplit[2]), 282, err_msg="Failed # marker sets test")
        np.testing.assert_almost_equal(float(lineSplit[9]), 99.64, decimal=2, err_msg="Failed completeness test")
        np.testing.assert_almost_equal(float(lineSplit[10]), 6.02, decimal=2, err_msg="Failed contamination test")
