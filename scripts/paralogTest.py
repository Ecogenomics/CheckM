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
import sys
import argparse
from collections import defaultdict
import time
import random

from lib.img import IMG
from Bio import Phylo

import numpy as np

class ParalogTest(object):
    def __init__(self):
        pass
    
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
            
            tree = Phylo.read(os.path.join(geneTreeDir, f), 'newick')
            
            taxa = tree.root.get_terminals()
            print '  Genes in tree: ' + str(len(taxa))
            
            # root tree with archaeal genomes
            print '  Attempting to rerooting tree between archaea and bacteria.'
            outgroup = []
            genomeIds = set()
            for t in taxa:
                try:
                    genomeId = t.name.split('_')[0]
                    #if '|' in genomeId:
                    #    genomeId = genomeId.split('|')[1]
                    genomeIds.add(genomeId)
                    domain =  metadata[genomeId]['taxonomy'][0].lower()
                except:
                    print t.name
                    sys.exit()
                    
                if domain == 'archaea':
                    outgroup.append(t)
                    
            print '  Genomes in tree: ' + str(len(genomeIds))
                   
            if not tree.is_monophyletic(outgroup):
                print '    Archaea are paraphyletic. Rerooting at midpoint.'
                tree.root_at_midpoint()
            else:
                print '    Archaea are monophyletic. Rerooting between archaea and bacteria.'
                tree.root_with_outgroup(outgroup)
                
            # get species name of each taxa
            leafNodeToSpeciesName = {}
            for t in taxa:
                genomeId = t.name.split('_')[0]
                genus = metadata[genomeId]['taxonomy'][5]
                sp = metadata[genomeId]['taxonomy'][6].lower()
                
                leafNodeToSpeciesName[t.name] = genus + ' ' + sp
                         
            # calculate intra-species phylogenetic distance
            if False:
                start = time.time()
                print '  Calculating intra-species distances.'
                
                intraSpeciesPairs = []
                for i in xrange(0, len(taxa)):
                    speciesI = leafNodeToSpeciesName[taxa[i].name]
                    for j in xrange(i+1, len(taxa)):
                        speciesJ = leafNodeToSpeciesName[taxa[j].name]
                        
                        if speciesI == speciesJ:
                            #intraSpeciesDist.append(tree.distance(taxa[i], taxa[j]))
                            intraSpeciesPairs.append([taxa[i], taxa[j]])
                            
                intraSpeciesDist = []
                random.shuffle(intraSpeciesPairs)
                for i in xrange(0, 100):
                    t1, t2 = intraSpeciesPairs[i]
                    intraSpeciesDist.append(tree.distance(t1, t2))
                            
                intraSpeciesDist = np.array(intraSpeciesDist) 
                
                intraSpeciesDist.sort()        
                percentile = intraSpeciesDist[int(0.75*len(intraSpeciesDist))]
                print '  Mean: %.3f, std: %.3f, median: %.5f, 80th percentile: %.4f' % (np.mean(intraSpeciesDist), np.std(intraSpeciesDist), np.median(intraSpeciesDist), percentile)
                end = time.time()
                print '  Time: ' + str(end - start)
  
            # find all paralogous genes
            print '  Finding paralogous genes.'
            
            paralogs = defaultdict(set)
            for i in xrange(0, len(taxa)):
                genomeId = taxa[i].name.split('_')[0]
                for j in xrange(i+1, len(taxa)):
                    # genes from the same genome paralogs, but we filter out
                    # those that are identical (distance of 0 on the tree) to 
                    # speed up computation and because these clearly do not 
                    # adversely effect phylogenetic inference
                    if genomeId == taxa[j].name.split('_')[0] and tree.distance(taxa[i], taxa[j]) > 0:
                        paralogs[genomeId].add(taxa[i])
                        paralogs[genomeId].add(taxa[j])
                        
            print '    Paralogous genes: ' + str(len(paralogs))
               
            # check if paralogous genes are conspecific
            print '  Determining if paralogous genes are conspecific.'
            nonConspecificGenomes = []
            for genomeId, taxa in paralogs.iteritems():
                lcaNode = tree.common_ancestor(taxa)
                
                children = lcaNode.get_terminals()
                species = set()
                for child in children:
                    genomeId = child.name.split('_')[0]
                    
                    genus = metadata[genomeId]['taxonomy'][5]
                    sp = metadata[genomeId]['taxonomy'][6].lower()
                    if sp != '' and sp != 'unclassified':
                        species.add(genus + ' ' + sp)
                        
                if len(species) > 1:
                    nonConspecificGenomes.append(genomeId)
            
            if len(nonConspecificGenomes) > acceptPer*len(genomeIds):
                filteredGeneTrees += 1
                print '  Tree is not conspecific for the following genome: ' + str(nonConspecificGenomes)
            else:
                retainedGeneTrees += 1
                
                if len(nonConspecificGenomes) > 1:
                    print '  An acceptable number of genomes are not conspecific: ' + str(nonConspecificGenomes)
                else:
                    print '  Tree is conspecific.'
                    
                os.system('cp ' + os.path.join(geneTreeDir, f) + os.path.join(outputDir, f))
                
                alignFile = f.replace('.tre', '.faa')
                os.system('cp ' + os.path.join(geneTreeDir, alignFile) + os.path.join(outputDir, alignFile))
                
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
