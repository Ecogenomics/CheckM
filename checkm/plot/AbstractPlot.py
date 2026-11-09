###############################################################################
#
# AbstractPlot.py - Abstract base class for plotting.
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

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.transforms as mtransforms

import matplotlib as mpl

import numpy as np


class AbstractPlot(FigureCanvas):
    '''
    Abstract base class for plotting.
    '''

    def __init__(self, options):
        self.options = options

        # Global plot settings
        mpl.rcParams['font.size'] = self.options.font_size
        mpl.rcParams['axes.titlesize'] = self.options.font_size
        mpl.rcParams['axes.labelsize'] = self.options.font_size
        mpl.rcParams['xtick.labelsize'] = self.options.font_size
        mpl.rcParams['ytick.labelsize'] = self.options.font_size
        mpl.rcParams['legend.fontsize'] = self.options.font_size
        mpl.rcParams['svg.fonttype'] = 'none'

        self.fig = Figure(facecolor='white', dpi=options.dpi)

        FigureCanvas.__init__(self, self.fig)

        self.axesColour = (0.5, 0.5, 0.5)

    def histogramBins(self, binWidth, maxValue=1.0):
        """Edges of bins with the given width that start at 0 and cover [0, maxValue]."""
        if binWidth <= 0:
            raise ValueError('Histogram bin width must be positive: %s' % str(binWidth))

        # small tolerance so floating point error does not create an extra bin
        numBins = max(int(np.ceil(maxValue / binWidth - 1e-6)), 1)
        return np.linspace(0, numBins * binWidth, numBins + 1)

    def percentWeights(self, data):
        """Histogram weights giving the percentage of data points in each bin."""
        return np.full(len(data), 100.0 / len(data))

    def savePlot(self, filename, dpi=300):
        imgFormat = filename[filename.rfind('.') + 1:len(filename)]
        if imgFormat in ['png', 'pdf', 'ps', 'eps', 'svg']:
            self.fig.savefig(filename, format=imgFormat, dpi=dpi,
                             facecolor='white', edgecolor='white', bbox_inches='tight')
        else:
            raise ValueError('Unsupported image format: %s' % imgFormat)

    def yLabelExtents(self, labels, fontSize, rotation=0):
        self.fig.clear()

        tempAxes = self.fig.add_axes([0, 0, 1.0, 1.0])
        tempAxes.set_yticks(np.arange(len(labels)))
        yLabels = tempAxes.set_yticklabels(
            labels, size=fontSize, rotation=rotation)

        bboxes = []
        for label in yLabels:
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.transformed(self.fig.transFigure.inverted())
            bboxes.append(bboxi)
        yLabelBounds = mtransforms.Bbox.union(bboxes)

        self.fig.clear()

        return yLabelBounds
