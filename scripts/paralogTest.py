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

from checkm.lib.img import IMG
import dendropy

from rerootTree import RerootTree

class ParalogTest(object):
    def __init__(self):
        pass
    
    def __patristicDist(self, tree, taxa1, taxa2):
        mrca = tree.mrca(taxon_labels=[taxa1.taxon.label, taxa2.taxon.label])
        
        if mrca.parent_node == None:
            # MRCA is the root of the tree
            return taxa1.distance_from_root() + taxa2.distance_from_root()
        else:
        
            dist = taxa1.edge_length
            parentNode = taxa1.parent_node
            while parentNode != mrca:                  
                dist += parentNode.edge_length    
                parentNode = parentNode.parent_node
    
                
            dist += taxa2.edge_length
            parentNode = taxa2.parent_node
            while parentNode != mrca:                  
                dist += parentNode.edge_length
                parentNode = parentNode.parent_node
                
            return dist

    def run(self, geneTreeDir, acceptPer, extension, outputDir):
        # make sure output directory is empty
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        files = os.listdir(outputDir)
        for f in files:
            os.remove(os.path.join(outputDir, f))

        img = IMG()
        metadata = img.genomeMetadata()

        files = os.listdir(geneTreeDir)
        print 'Identifying gene trees with only conspecific paralogous genes:'
        filteredGeneTrees = 0
        retainedGeneTrees = 0
        for f in files:
            if not f.endswith(extension):
                continue

            geneId = f[0:f.find('.')]
            print '  Testing gene tree: ' + geneId

            tree = dendropy.Tree.get_from_path(os.path.join(geneTreeDir, f), schema='newick', as_rooted=False, preserve_underscores=True)

            taxa = tree.leaf_nodes()
            numTaxa = len(taxa)
            print '  Genes in tree: ' + str(numTaxa)

            # root tree with archaeal genomes
            rerootTree = RerootTree()
            rerootTree.reroot(tree)
            
            # get species name of each taxa
            leafNodeToSpeciesName = {}
            for t in taxa:
                genomeId = t.taxon.label.split('|')[0]
                genus = metadata[genomeId]['taxonomy'][5]
                sp = metadata[genomeId]['taxonomy'][6].lower()

                leafNodeToSpeciesName[t.taxon.label] = genus + ' ' + sp
                
            # find all paralogous genes
            print '  Finding paralogous genes.'

            paralogs = defaultdict(set)
            for i in xrange(0, len(taxa)):
                genomeId = taxa[i].taxon.label.split('|')[0]
                for j in xrange(i+1, len(taxa)):
                    # genes from the same genome are paralogs, but we filter out
                    # those that are identical (distance of 0 on the tree) to
                    # speed up computation and because these clearly do not
                    # adversely effect phylogenetic inference
                    if genomeId == taxa[j].taxon.label.split('|')[0] and self.__patristicDist(tree, taxa[i], taxa[j]) > 0:
                        paralogs[genomeId].add(taxa[i].taxon.label)
                        paralogs[genomeId].add(taxa[j].taxon.label)
                        
            print '    Paralogous genes: ' + str(len(paralogs))

            # check if paralogous genes are conspecific
            print '  Determining if paralogous genes are conspecific.'
            nonConspecificGenomes = []
            for genomeId, taxaLabels in paralogs.iteritems():
                lcaNode = tree.mrca(taxon_labels = taxaLabels)

                children = lcaNode.leaf_nodes()
                species = set()
                for child in children:
                    childGenomeId = child.taxon.label.split('|')[0]

                    genus = metadata[childGenomeId]['taxonomy'][5]
                    sp = metadata[childGenomeId]['taxonomy'][6].lower()
                    if sp != '' and sp != 'unclassified':
                        species.add(genus + ' ' + sp)

                if len(species) > 1:
                    nonConspecificGenomes.append(genomeId)

            if len(nonConspecificGenomes) > acceptPer*numTaxa:
                filteredGeneTrees += 1
                print '  Tree is not conspecific for the following genome: ' + str(nonConspecificGenomes)
            else:
                retainedGeneTrees += 1

                if len(nonConspecificGenomes) > 1:
                    print '  An acceptable number of genomes are not conspecific: ' + str(nonConspecificGenomes)
                else:
                    print '  Tree is conspecific.'

                os.system('cp ' + os.path.join(geneTreeDir, f) + ' ' + os.path.join(outputDir, f))

            print ''

        print 'Filtered gene trees: ' + str(filteredGeneTrees)
        print 'Retained gene trees: ' + str(retainedGeneTrees)

if __name__ == '__main__':
    print 'ParalogTest v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gene_tree_dir', help='directory containing gene trees to test')
    parser.add_argument('-a', '--accept_per', help='percentage of non-conspecific genomes allowed', type=float, default=0.003)
    parser.add_argument('-x', '--extension', help='extension of tree files to process', default = '.tre')
    parser.add_argument('-o', '--output_dir', help='output directory for retained gene trees', default = './data/gene_trees_conspecific/')

    args = parser.parse_args()

    paralogTest = ParalogTest()
    paralogTest.run(args.gene_tree_dir, args.accept_per, args.extension, args.output_dir)
