###############################################################################
#
# nxPlot.py - Create a Nx-plot.
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

import sys
import logging

import numpy as np

from matplotlib.ticker import FuncFormatter

from checkm.plot.AbstractPlot import AbstractPlot

from checkm.util.seqUtils import readFasta


def nxThresholds(stepSize):
    """Nx thresholds from 0 to 1 (inclusive) in increments of stepSize."""
    x = np.arange(0, 1.0, stepSize)
    x = x[x < 1.0 - 1e-6 * stepSize]  # avoid near duplicate of 1.0 due to floating point error
    return np.append(x, 1.0)


def kbpLabel(seqLen, _pos=None):
    """Format a sequence length in bp as a label in kbp."""
    label = '%.1f' % (float(seqLen) / 1000)
    if label.endswith('.0'):
        label = label[0:-2]  # remove trailing zero

    return label


class NxPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def calculateNx(self, x, seqs):
        """Calculate Nx for each threshold in x (fractions of the total sequence length).

        Nx is the length of the shortest sequence such that sequences
        of this length or longer contain at least x of all bases.
        """
        seqLens = np.array(sorted([len(seq) for seq in seqs.values()], reverse=True))
        if len(seqLens) == 0 or seqLens.sum() == 0:
            return [0] * len(x)

        cumSum = np.cumsum(seqLens)
        indices = np.searchsorted(cumSum, np.asarray(x) * cumSum[-1], side='left')
        indices = np.minimum(indices, len(seqLens) - 1)

        return [int(v) for v in seqLens[indices]]

    def plot(self, fastaFile):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axes = self.fig.add_subplot(111)

        # calculate Nx
        if not 0 < self.options.step_size <= 1:
            logging.getLogger('timestamp').error('The step size must be > 0 and <= 1: %s' % str(self.options.step_size))
            sys.exit(1)

        seqs = readFasta(fastaFile)
        x = nxThresholds(self.options.step_size)
        nx = self.calculateNx(x, seqs)

        # Create plot
        axes.plot(x, nx, 'ko-', ms=4)
        axes.set_xlabel('Nx')
        axes.set_ylabel('Sequence length (kbp)')

        # Change sequence lengths from bp to kbp
        axes.yaxis.set_major_formatter(FuncFormatter(kbpLabel))

        # Prettify plot
        for line in axes.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axes.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in axes.spines.items():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)

        self.fig.tight_layout(pad=1)
        self.draw()
