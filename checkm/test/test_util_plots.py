###############################################################################
#
# test_util_plots.py - unit tests for plotting helpers
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
import tempfile
import unittest
import warnings
from argparse import Namespace

import numpy as np

import matplotlib
matplotlib.use('Agg')

from checkm.plot.nxPlot import NxPlot, nxThresholds, kbpLabel
from checkm.plot.lengthHistogram import LengthHistogram
from checkm.util.seqUtils import writeFasta, calculateN50


def _options(**kwargs):
    options = Namespace(font_size=8, dpi=50, width=4, height=3, step_size=0.05)
    for k, v in kwargs.items():
        setattr(options, k, v)
    return options


class VerifyNxPlot(unittest.TestCase):
    def setUp(self):
        rng = random.Random(0)
        self.seqs = {'s%d' % i: 'A' * rng.randint(100, 50000) for i in range(57)}
        self.plot = NxPlot(_options())

    def testThresholds(self):
        for stepSize in np.linspace(0.001, 0.5, 500):
            x = nxThresholds(stepSize)
            self.assertEqual(x[0], 0)
            self.assertEqual(x[-1], 1.0)
            self.assertTrue(np.all(np.diff(x) > 0))
            self.assertTrue(np.all(np.diff(x)[:-1] <= stepSize * (1 + 1e-9)))

        self.assertEqual(len(nxThresholds(0.1)), 11)
        self.assertEqual(len(nxThresholds(0.05)), 21)
        self.assertEqual(len(nxThresholds(0.3)), 5)  # 0, 0.3, 0.6, 0.9, 1.0
        self.assertEqual(list(nxThresholds(1.0)), [0, 1.0])

    def testCalculateNx(self):
        seqLens = sorted([len(s) for s in self.seqs.values()], reverse=True)

        nx = self.plot.calculateNx([0, 0.5, 1.0], self.seqs)
        self.assertEqual(nx[0], seqLens[0])
        self.assertEqual(nx[1], calculateN50(seqLens))
        self.assertEqual(nx[2], seqLens[-1])

        # brute force check
        x = nxThresholds(0.013)
        nx = self.plot.calculateNx(x, self.seqs)
        self.assertEqual(len(nx), len(x))
        total = sum(seqLens)
        for xi, nxi in zip(x, nx):
            expected = None
            cumSum = 0
            for seqLen in seqLens:
                cumSum += seqLen
                if cumSum >= xi * total:
                    expected = seqLen
                    break
            self.assertEqual(nxi, expected)

        self.assertEqual(self.plot.calculateNx([0, 0.5, 1.0], {}), [0, 0, 0])

    def testPlotStepSizes(self):
        """Nx plots can be created for any step size."""
        tmpDir = tempfile.mkdtemp()
        try:
            fastaFile = os.path.join(tmpDir, 'bin.fna')
            writeFasta(self.seqs, fastaFile)

            for stepSize in [0.001, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.13, 0.25, 0.3, 0.33, 0.5, 1.0]:
                plot = NxPlot(_options(step_size=stepSize))
                with warnings.catch_warnings():
                    warnings.simplefilter('error', UserWarning)
                    plot.plot(fastaFile)
                plot.savePlot(os.path.join(tmpDir, 'nx.png'), dpi=50)
        finally:
            shutil.rmtree(tmpDir)

    def testKbpLabel(self):
        self.assertEqual(kbpLabel(0), '0')
        self.assertEqual(kbpLabel(10000), '10')
        self.assertEqual(kbpLabel(10500), '10.5')
        self.assertEqual(kbpLabel(100000), '100')


class VerifyHistograms(unittest.TestCase):
    def setUp(self):
        self.plot = LengthHistogram(_options())

    def testHistogramBins(self):
        """Bins cover the full range [0, 1] including the top bin."""
        for binWidth in [0.01, 0.02, 0.025, 0.03, 0.05, 0.07, 0.1, 0.3, 1.0]:
            bins = self.plot.histogramBins(binWidth)
            self.assertEqual(bins[0], 0)
            self.assertGreaterEqual(bins[-1], 1.0)
            self.assertLess(bins[-1] - binWidth, 1.0)
            self.assertTrue(np.allclose(np.diff(bins), binWidth))

        bins = self.plot.histogramBins(0.01)
        self.assertEqual(len(bins), 101)
        self.assertEqual(bins[-1], 1.0)

        # a value in the top bin is counted
        counts, _ = np.histogram([0.995, 1.0], bins=bins)
        self.assertEqual(counts[-1], 2)

        with self.assertRaises(ValueError):
            self.plot.histogramBins(0)

    def testPercentWeights(self):
        data = [0.1, 0.2, 0.2, 0.9]
        for binWidth in [0.01, 0.05, 0.1]:
            counts, _ = np.histogram(data, bins=self.plot.histogramBins(binWidth),
                                     weights=self.plot.percentWeights(data))
            self.assertAlmostEqual(counts.sum(), 100.0)
            self.assertAlmostEqual(counts.max(), 50.0)


class VerifySavePlot(unittest.TestCase):
    def testUnknownFormat(self):
        plot = LengthHistogram(_options())
        with self.assertRaises(ValueError):
            plot.savePlot('plot.unknown')


if __name__ == "__main__":
    unittest.main()
