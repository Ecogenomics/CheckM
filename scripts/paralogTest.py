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

__prog_desc__ = 'determine if paralogous genes are conspecific'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import argparse
from collections import defaultdict

from lib.img import IMG
from Bio import Phylo

class ParalogTest(object):
    def __init__(self):
        pass
    
    def treeLCA(self, tree, taxa):
        """ Find lowest common ancestor of a set of taxa."""
        lca = tree.root
        for node in tree.trace(list(taxa)[0], tree.root):
            bLCA = True
            for l in taxa:
                if l not in node.get_terminals():
                    bLCA = False
                    break
    
            if bLCA:
                lca = node
    
        return lca
    
    def run(self, geneTreeDir, metadataFile, extension):
        metadata = img.genomeMetadata()
        
        files = os.listdir(geneTreeDir)
        print 'Genes with only conspecific paralogous genes:'
        for f in files:
            if not f.endswith(extension):
                continue
            
            tree = Phylo.read(os.path.join(geneTreeDir, f), 'newick')
              
            # find all paralogous genes
            taxa = tree.root.get_terminals()
            
            paralogs = defaultdict(set)
            for i in xrange(0, len(taxa)):
                genomeId = taxa[i].name.split('_')[0]
                for j in xrange(i+1, len(taxa)):
                    if genomeId == taxa[j].name.split('_')[0]:
                        paralogs[genomeId].add(taxa[j])
               
            # check if paralogous genes are conspecific
            for genomeId, taxa in paralogs.iteritems():
                lcaNode = self.treeLCA(tree, taxa)
                
                children = lcaNode.get_terminals()
                species = set()
                for child in children:
                    genomeId = child[0].name.split('_')[0]
                    
                    sp = metadata[genomeId]['taxonomy']['species'].lower()
                    if sp != '' and sp != 'unclassified':
                        species.add(sp)
                        
            # sanity checking
            if len(species) > 1:
                print species

            # report status of each gene
            geneId = f[0:f.rfind('.')]
            if len(species) == 1:
                print geneId
            
if __name__ == '__main__':
    print 'ParalogTest v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'
  
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gene_tree_dir', help='directory containing gene trees to test')
    parser.add_argument('-x', '--extension', help='extension of tree files to process', default = '.tre')
    
    args = parser.parse_args()
    
    paralogTest = ParalogTest()
    paralogTest.run(args.gene_tree_dir, args.extension)
