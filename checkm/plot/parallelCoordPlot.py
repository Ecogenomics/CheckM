###############################################################################
#
# parallelCoordPlot.py - Create a parallel coordinate plot. 
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

import numpy as np

from AbstractPlot import AbstractPlot

class parallelCoordPlot(AbstractPlot):
	def __init__(self, options):
		AbstractPlot.__init__(self, options)

	def plot(self, seqStats, converageStats):       
		# Set size of figure
		self.fig.clear()
		self.fig.set_size_inches(self.options.width, self.options.height)
		#heightBottomLabels = self.options.fig_padding   # inches
		#widthSideLabel = self.options.fig_padding       # inches 
		#axes = self.fig.add_axes([widthSideLabel/self.options.width,heightBottomLabels/self.options.height,\
		#																1.0-(widthSideLabel+self.options.fig_padding )/self.options.width,\
		#																1.0-(heightBottomLabels+self.options.fig_padding )/self.options.height])
																		
		# create data points for each sequence
		data = []
		for seqId, seqStat in seqStats:
			data.append([seqStat['GC']] + coverageStats[seqId].values())
			
		dims = len(data[0])
		x = range(dims)
		fig, axes = self.fig.add_subplots(1, dims-1, sharey=False)

		if style is None:
			style = ['r-']*len(data)

		# Calculate the limits on the data
		min_max_range = list()
		for m in zip(*data):
			mn = min(m)
			mx = max(m)
			if mn == mx:
				mn -= 0.5
				mx = mn + 1.
			r  = float(mx - mn)
			min_max_range.append((mn, mx, r))

		# Normalize the data sets
		norm_data_sets = list()
		for ds in data:
			nds = [(value - min_max_range[dimension][0]) / 
					min_max_range[dimension][2] 
					for dimension,value in enumerate(ds)]
			norm_data_sets.append(nds)
		data = norm_data_sets

		# Plot the datasets on all the subplots
		for i, ax in enumerate(axes):
			for dsi, d in enumerate(data):
				ax.plot(x, d, style[dsi])
			ax.set_xlim([x[i], x[i+1]])

		# Set the x axis ticks 
		for dimension, (axx,xx) in enumerate(zip(axes, x[:-1])):
			axx.xaxis.set_major_locator(ticker.FixedLocator([xx]))
			ticks = len(axx.get_yticklabels())
			labels = list()
			step = min_max_range[dimension][2] / (ticks - 1)
			mn   = min_max_range[dimension][0]
			for i in xrange(ticks):
				v = mn + i*step
				labels.append('%4.2f' % v)
			axx.set_yticklabels(labels)

		# Move the final axis' ticks to the right-hand side
		axx = plt.twinx(axes[-1])
		dimension += 1
		axx.xaxis.set_major_locator(ticker.FixedLocator([x[-2], x[-1]]))
		ticks = len(axx.get_yticklabels())
		step = min_max_range[dimension][2] / (ticks - 1)
		mn   = min_max_range[dimension][0]
		labels = ['%4.2f' % (mn + i*step) for i in xrange(ticks)]
		axx.set_yticklabels(labels)

		# Stack the subplots 
		plt.subplots_adjust(wspace=0)
	
		# Prettify plot     
		for a in axes.yaxis.majorTicks:
			a.tick1On=True
			a.tick2On=False
				
		for a in axes.xaxis.majorTicks:
			a.tick1On=True
			a.tick2On=False
			
		for line in axes.yaxis.get_ticklines(): 
			line.set_color(self.axesColour)
				
		for line in axes.xaxis.get_ticklines(): 
			line.set_color(self.axesColour)
			
		for loc, spine in axes.spines.iteritems():
			if loc in ['right','top']:
				spine.set_color('none') 
			else:
				spine.set_color(self.axesColour)
						  
		self.draw()