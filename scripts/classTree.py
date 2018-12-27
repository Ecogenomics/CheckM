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

__prog_desc__ = 'create a reference genome tree suitable for collapsing'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import random
from collections import defaultdict

import dendropy

from checkm.lib.img import IMG

class ClassTree(object):
    def __init__(self):
        self.img = IMG()
    
    def run(self):
        tree = dendropy.Tree.get_from_path('../data/genome_tree/genome_tree_prok.refpkg/genome_tree.final.tre', schema='newick', as_rooted=True, preserve_underscores=True)
        
        metadata = self.img.genomeMetadata()
        
        # relabel taxa
        for leaf in tree.leaf_nodes():
            genomeId = leaf.taxon.label.replace('IMG_', '')
            classT = metadata[genomeId]['taxonomy'][2]
            newLeafLabel = classT + '_' + genomeId
            leaf.taxon.label = newLeafLabel
            
        # relabel internal nodes
        for node in tree.internal_nodes():
            uid, taxaStr, bootstrap = node.label.split('|')
            
            if bootstrap:
                node.label = uid + ':' + str(node.edge_length) + '[' + str(int(float(bootstrap)*100 + 0.5)) + ']'
            else:
                node.label = uid + ':' + str(node.edge_length)

        # write out reduced tree
        tree.write_to_path('./experiments/classTree.tre', schema='newick', suppress_rooting=True, suppress_edge_lengths=True, unquoted_underscores=True)
        tree.write_to_path('./experiments/classTree_no_internal.tre', schema='newick', suppress_rooting=True, suppress_edge_lengths=True, unquoted_underscores=True, suppress_internal_node_labels=True)

if __name__ == '__main__':
    print('ClassTree v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    classTree = ClassTree()
    classTree.run()
