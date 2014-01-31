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

__prog_desc__ = 'propogate internal labels from one tree to another'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse

import dendropy

class PropogateTaxonomicLabels(object):
    def __init__(self):
        pass

    def propogate(self, labelTree, derepFile, inputTree):
        # find primary taxa
        primaryTaxa = set()
        for line in open(derepFile):
            primaryTaxa.add(line.split()[0])
        
        # find internal nodes with labels
        for node in labelTree.internal_nodes():
            if node.label != None: 
                labels = []
                for leaf in node.leaf_nodes():
                    if leaf.taxon.label in primaryTaxa:
                        labels.append(leaf.taxon.label)
                
                if len(labels) > 1:
                    # find MRCA of leaves and amend label
                    taxaStr = node.label.split(':')[0].replace(' ', '')
                    mrca = inputTree.mrca(taxon_labels = labels)
                    mrca.label = taxaStr + '|' + mrca.label
                else:
                    # label is undefined in new tree as it belongs to a single genome
                    pass

    def run(self, labelTreeFile, derepFile, inputTreeFile):
        labelTree = dendropy.Tree.get_from_path(labelTreeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        inputTree = dendropy.Tree.get_from_path(inputTreeFile, schema='newick', as_rooted=True, preserve_underscores=True)

        self.propogate(labelTree, derepFile, inputTree)

        inputTree.write_to_path(inputTreeFile, schema='newick', suppress_rooting=True, unquoted_underscores=True)

if __name__ == '__main__':
    print 'PropogateTaxonomicLabels v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('label_tree', help='tree with internal labels to propogate')
    parser.add_argument('derep_file', help='file indicating dereplicated nodes in input tree')
    parser.add_argument('input_tree', help='tree to decorate with new labels')

    args = parser.parse_args()

    propogateTaxonomicLabels = PropogateTaxonomicLabels()
    propogateTaxonomicLabels.run(args.label_tree, args.derep_file, args.input_tree)
