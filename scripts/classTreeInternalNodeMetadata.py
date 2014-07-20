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
Get marker set metadata for class tree.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

from checkm.treeParser import TreeParser
from checkm.markerSets import MarkerSet

class ClassTreeMetadata(object):
    def __init__(self):
        pass

    def run(self):
        # read internal nodes file
        metadata = {}
        for line in open('./experiments/classTree.internal_nodes.tsv'):
            uid, label = [x.strip() for x in line.split('\t')]
            metadata[uid] = label
        

        # read all lineage-specific marker genes
        treeParser = TreeParser()
        uniqueIdToLineageStatistics = treeParser.readNodeMetadata()
        for uid in metadata:
            stats = uniqueIdToLineageStatistics[uid]
            markerSet = MarkerSet(uid, 'NA', int(stats['# genomes']), eval(stats['marker set']))
            
            metadata[uid] += ' [%d, %d, %d]' % (stats['# genomes'], markerSet.numMarkers(), markerSet.numSets())
            
        # write out results
        fout = open('./experiments/classTree.internal_nodes.metadata.tsv', 'w')
        for uid, label in metadata.iteritems():
            fout.write(uid + '\t' + label + '\n')
        fout.close()
            
if __name__ == '__main__':
    classTreeMetadata = ClassTreeMetadata()
    classTreeMetadata.run()
