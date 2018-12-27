#!/usr/bin/env python

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

"""
Gather all phylogenetically informative HMMs into a single model file.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import argparse

class GetPhylogeneticHMMs(object):
    def __init__(self):
        pass

    def run(self, phyloHMMs, geneTrees, outputFile):
        # get list of phylogenetically informative marker genes
        files = os.listdir(geneTrees)
        
        markerIds = []
        for f in files:
            if f.endswith('.tre'):
                markerId = f[0:f.find('.')]
                markerIds.append(markerId)
                
        print('Identified %d phylogenetically informative marker genes.' % (len(markerIds)))
                
        # place all phylogenetically informative marker genes into a single model file
        fout = open(outputFile, 'w')
        for markerId in markerIds:
            for line in open(os.path.join(phyloHMMs, markerId + '.hmm')):
                fout.write(line)
        
        fout.close()
        
        print('HMMs written to: ' + outputFile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Gather all phylogenetically informative HMMs into a single model file.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('phylo_hmms', help='directory containing all potentially phylogenetically informative HMMs')
    parser.add_argument('gene_trees', help='directory containing phylogenetically informative gene trees')
    parser.add_argument('output_file', help='output HMM model file')

    args = parser.parse_args()

    getPhylogeneticHMMs = GetPhylogeneticHMMs()
    getPhylogeneticHMMs.run(args.phylo_hmms, args.gene_trees, args.output_file)
