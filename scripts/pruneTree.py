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

__prog_desc__ = 'prune taxa with identical sequences from tree'

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

class PruneTree(object):
    def __init__(self):
        pass
    
    def __readDuplicateTaxa(self, dupSeqFile):
        dupTaxa = []
        
        for line in open(dupSeqFile):
            lineSplit = line.split()
            
            for i in range(1, len(lineSplit)):
                dupTaxa.append(lineSplit[i].strip())
                
        return dupTaxa
                

    def run(self, dupSeqFile, inputTree, outputTree):
        # get list of taxa with duplicate sequences 
        dupTaxa = self.__readDuplicateTaxa(dupSeqFile)
        print('Pruing %d taxa.' % len(dupTaxa))
         
        # prune duplicate taxa from tree
        tree = dendropy.Tree.get_from_path(inputTree, schema='newick', as_rooted=True, preserve_underscores=True)

        tree.prune_taxa_with_labels(dupTaxa)

        tree.write_to_path(outputTree, schema='newick', suppress_rooting=True)

if __name__ == '__main__':
    print('RerootTree v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('duplicate_seq_file', help='file indicating deplicate sequences as determine with seqmagick')
    parser.add_argument('input_tree', help='tree to prunt')
    parser.add_argument('output_tree', help='output tree')

    args = parser.parse_args()

    PruneTree = PruneTree()
    PruneTree.run(args.duplicate_seq_file, args.input_tree, args.output_tree)
