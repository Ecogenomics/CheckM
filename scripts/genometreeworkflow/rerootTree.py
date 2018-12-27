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

__prog_desc__ = 'attempts to reroot a tree between bacteria and archaea'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse

from checkm.util.img import IMG
import dendropy

class RerootTree(object):
    def __init__(self):
        img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        self.metadata = img.genomeMetadata()

    def __genomeId(self, taxaLabel):
        if taxaLabel.startswith('IMG_'):
            genomeId = taxaLabel[4:]
        elif '|' in taxaLabel:
            # assume taxa are labeled with the format: genomeId_geneId
            genomeId = taxaLabel.split('|')[0]
        else:
            genomeId = taxaLabel

        return genomeId
    
    def __isMonophyletic(self, tree, group):
        tree.update_splits()
        
        groupSplit = tree.taxon_set.get_taxa_bitmask(labels=group)
        return groupSplit in tree.split_edges

    def reroot(self, tree):
        taxa = tree.leaf_nodes()

        # root tree with archaeal genomes
        print('  Attempting to rerooting tree between archaea and bacteria.')
        outgroup = []
        genomeIds = set()
        for t in taxa:
            try:
                genomeId = self.__genomeId(t.taxon.label)
                genomeIds.add(genomeId)
                domain = self.metadata[genomeId]['taxonomy'][0].lower()
            except:
                print('[Error] Missing IMG metadata for: ' + t.taxon.label)
                sys.exit()

            if domain == 'archaea':
                outgroup.append(t.taxon.label)

        if outgroup == [] or len(outgroup) == len(taxa) or not self.__isMonophyletic(tree, outgroup):
            print('    Archaea are paraphyletic or tree is exclusively bacteria or archaea. Rerooting at midpoint.')
            tree.reroot_at_midpoint(update_splits=True)
        else:
            print('    Archaea are monophyletic. Rerooting between archaea and bacteria.')
            mrca = tree.mrca(taxon_labels=outgroup)
            tree.reroot_at_edge(mrca.edge, length1 = 0.5*mrca.edge_length, length2 = 0.5*mrca.edge_length, update_splits=True)

    def run(self, inputTree, outputTree):
        tree = dendropy.Tree.get_from_path(inputTree, schema='newick', as_rooted=False, preserve_underscores=True)

        self.reroot(tree)

        tree.write_to_path(outputTree, schema='newick', suppress_rooting=True, unquoted_underscores=True)

if __name__ == '__main__':
    print('RerootTree v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_tree', help='tree to reroot')
    parser.add_argument('output_tree', help='rerooted tree')

    args = parser.parse_args()

    rerootTree = RerootTree()
    rerootTree.run(args.input_tree, args.output_tree)
