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

__prog_desc__ = 'complete workflow for inferring a genome tree for CheckM'

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

from phylogeneticInferenceGenes import PhylogeneticInferenceGenes
from makeTrees import MakeTrees
from paralogTest import ParalogTest
from consistencyTest import ConsistencyTest
from getPhylogeneticHMMs import GetPhylogeneticHMMs
from inferGenomeTree import InferGenomeTree
from rerootTree import RerootTree
from bootstrapTree import BootstrapTree
from propogateTaxonomicLabels import PropogateTaxonomicLabels

class GenomeTreeWorkflow(object):
    def __init__(self, outputDir):
        if os.path.exists(outputDir):
            print '[Error] Output directory already exists: ' + outputDir
            sys.exit(0)
        else:
            os.makedirs(outputDir)
            
        self.__checkForHMMER()
        self.__checkForFastTree()
        self.__checkForSeqMagick()
        self.__checkForTax2Tree()
            
        self.hmmDir = os.path.join(outputDir, 'phylo_hmms')
        self.alignmentDir = os.path.join(outputDir, 'gene_alignments')
        self.geneTreeDir = os.path.join(outputDir, 'gene_trees')
        self.conspecificGeneTreeDir = os.path.join(outputDir, 'gene_trees_conspecific')
        self.finalGeneTreeDir = os.path.join(outputDir, 'gene_trees_final')
        self.bootstrapDir = os.path.join(outputDir, 'bootstrap')
        
        self.consistencyOut = os.path.join(outputDir, 'genome_tree.consistency.tsv')
        self.concatenatedAlignFile = os.path.join(outputDir, 'genome_tree.concatenated.faa')
        self.derepConcatenatedAlignFile = os.path.join(outputDir, 'genome_tree.concatenated.derep.fasta')
        self.genomeTreeOut = os.path.join(outputDir, 'genome_tree.tre')
        self.genomeTreeRootedOut = os.path.join(outputDir, 'genome_tree.rooted.tre')
        self.genomeTreeDecoratedOut = os.path.join(outputDir, 'genome_tree.decorated.tre')
        self.genomeTreeDerepOut = os.path.join(outputDir, 'genome_tree.derep.tre')
        self.genomeTreeDerepRootedOut = os.path.join(outputDir, 'genome_tree.derep.rooted.tre')
        self.genomeTreeBootstrapOut = os.path.join(outputDir, 'genome_tree.bs.tre')
        self.genomeTreeFinalOut = os.path.join(outputDir, 'genome_tree.final.tre')
        self.genomeTreeTaxonomyOut = os.path.join(outputDir, 'genome_tree.taxonomy.tsv')
        self.phyloHMMsOut = os.path.join(outputDir, 'phylo.hmm')
        self.derepSeqFile = os.path.join(outputDir, 'genome_tree.derep.txt')
                
        self.phyloUbiquity = 0.9
        self.phyloSingleCopy = 0.9
        self.paralogAcceptPer = 0.003
        self.consistencyAcceptPer = 0.86
        self.consistencyMinTaxa = 20
        
        # create output directories
        os.makedirs(self.hmmDir)
        os.makedirs(self.alignmentDir)
        os.makedirs(self.geneTreeDir)
        os.makedirs(self.conspecificGeneTreeDir)
        os.makedirs(self.finalGeneTreeDir)
        os.makedirs(self.bootstrapDir)
    
    def __checkForHMMER(self):
        """Check to see if HMMER is on the system path."""

        try:
            exit_status = os.system('hmmfetch -h > /dev/null')
        except:
            print "Unexpected error!", sys.exc_info()[0]
            raise
    
        if exit_status != 0:
            print "[Error] hmmfetch is not on the system path"
            sys.exit()
            
    def __checkForFastTree(self):
        """Check to see if FastTree is on the system path."""

        try:
            exit_status = os.system('FastTree 2> /dev/null')
        except:
            print "Unexpected error!", sys.exc_info()[0]
            raise
    
        if exit_status != 0:
            print "[Error] FastTree is not on the system path"
            sys.exit()
            
    def __checkForSeqMagick(self):
        """Check to see if seqmagick is on the system path."""

        try:
            exit_status = os.system('seqmagick -h > /dev/null')
        except:
            print "Unexpected error!", sys.exc_info()[0]
            raise
    
        if exit_status != 0:
            print "[Error] seqmagick is not on the system path"
            sys.exit()
            
    def __checkForTax2Tree(self):
        """Check to see if tax2tree is on the system path."""

        try:
            exit_status = os.system('nlevel -h > /dev/null')
        except:
            print "Unexpected error!", sys.exc_info()[0]
            raise
    
        if exit_status != 0:
            print "[Error] nlevel is not on the system path"
            sys.exit()

    def run(self, numThreads):  
        # identify genes suitable for phylogenetic inference
        print '--- Identifying genes suitable for phylogenetic inference ---'
        phylogeneticInferenceGenes = PhylogeneticInferenceGenes()
        phylogeneticInferenceGenes.run(self.phyloUbiquity, self.phyloSingleCopy, numThreads, self.alignmentDir, self.hmmDir)
    
        # infer gene trees
        print ''
        print '--- Inferring gene trees ---'
        makeTrees = MakeTrees()
        makeTrees.run(self.alignmentDir, self.geneTreeDir, '.aln.masked.faa', numThreads)
        
        # test gene trees for paralogs
        print ''
        print '--- Testing for paralogs in gene trees ---'
        paralogTest = ParalogTest()
        paralogTest.run(self.geneTreeDir, self.paralogAcceptPer, '.tre', self.conspecificGeneTreeDir)
        
        # test gene trees for consistency with IMG taxonomy
        print ''
        print '--- Testing taxonomic consistency of gene trees ---'
        consistencyTest = ConsistencyTest()
        consistencyTest.run(self.conspecificGeneTreeDir, '.tre', self.consistencyAcceptPer, self.consistencyMinTaxa, self.consistencyOut, self.finalGeneTreeDir)
        
        # gather phylogenetically informative HMMs into a single model file
        print ''
        print '--- Gathering phylogenetically informative HMMs ---'
        getPhylogeneticHMMs = GetPhylogeneticHMMs()
        getPhylogeneticHMMs.run(self.hmmDir, self.finalGeneTreeDir, self.phyloHMMsOut)
        
        # infer genome tree
        print ''
        print '--- Inferring full genome tree ---'
        inferGenomeTree = InferGenomeTree()
        inferGenomeTree.run(self.finalGeneTreeDir, self.alignmentDir, '.aln.masked.faa', self.concatenatedAlignFile, self.genomeTreeOut, self.genomeTreeTaxonomyOut)
        
        # root genome tree between archaea and bacteria
        print ''
        print '--- Rooting full genome tree ---'
        rerootTree = RerootTree()
        rerootTree.run(self.genomeTreeOut, self.genomeTreeRootedOut)
            
        # decorate genome tree with taxonomy using nlevel from tax2tree
        print ''
        print '--- Decorating full genome tree with taxonomic information using tax2tree ---'
        os.system('nlevel -t %s -m %s -o %s' % (self.genomeTreeRootedOut, self.genomeTreeTaxonomyOut, self.genomeTreeDecoratedOut))
        
        # dereplicate identical sequences   
        print ''
        print '--- Identifying duplicate sequences ---'
        os.system('seqmagick convert --deduplicate-sequences --deduplicated-sequences-file ' + self.derepSeqFile + ' ' + self.concatenatedAlignFile + ' ' + self.derepConcatenatedAlignFile)
        
        # infer dereplicated genome tree 
        print ''
        print '--- Inferring dereplicated genome tree ---'
        outputLog = self.genomeTreeDerepOut[0:self.genomeTreeDerepOut.rfind('.')] + '.log'
        cmd = 'FastTreeMP -nosupport -wag -gamma -log ' + outputLog + ' ' + self.derepConcatenatedAlignFile + ' > ' + self.genomeTreeDerepOut
        os.system(cmd)
        
        # root genome tree between archaea and bacteria
        print ''
        print '--- Rooting dereplicated genome tree ---'
        rerootTree = RerootTree()
        rerootTree.run(self.genomeTreeDerepOut, self.genomeTreeDerepRootedOut)
                
        # calculate bootstraps for genome tree   
        print ''
        print '--- Calculating bootstrap support ---'
        bootstrapTree = BootstrapTree()
        bootstrapTree.run(self.bootstrapDir, self.genomeTreeDerepRootedOut, self.concatenatedAlignFile, 100, numThreads, self.genomeTreeBootstrapOut)
        
        os.system('cp ' + self.genomeTreeBootstrapOut + ' ' + self.genomeTreeFinalOut)
        
        # propogate taxonomic labels onto dereplicated tree 
        print ''
        print '--- Propagating taxonomic labels to dereplicated tree ---'
        propogateTaxonomicLabels = PropogateTaxonomicLabels()
        propogateTaxonomicLabels.run(self.genomeTreeDecoratedOut, self.derepSeqFile, self.genomeTreeFinalOut)
        
if __name__ == '__main__':
    print 'GenomeTreeWorkflow v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)

    args = parser.parse_args()

    genomeTreeWorkflow = GenomeTreeWorkflow(args.output_dir)
    genomeTreeWorkflow.run(args.threads)