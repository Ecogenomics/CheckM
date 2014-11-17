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

from genometreeworkflow.phylogeneticInferenceGenes import PhylogeneticInferenceGenes
from genometreeworkflow.makeTrees import MakeTrees
from genometreeworkflow.paralogTest import ParalogTest
from genometreeworkflow.consistencyTest import ConsistencyTest
from genometreeworkflow.getPhylogeneticHMMs import GetPhylogeneticHMMs
from genometreeworkflow.inferGenomeTree import InferGenomeTree
from genometreeworkflow.rerootTree import RerootTree
from genometreeworkflow.decorateTree import DecorateTree

class GenomeTreeWorkflow(object):
    def __init__(self, outputDir):
        if False:
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
        self.treeOut = os.path.join(outputDir, 'genome_tree.tre')
        self.treeRootedOut = os.path.join(outputDir, 'genome_tree.rooted.tre')
        self.treeTaxonomyOut = os.path.join(outputDir, 'genome_tree.taxonomy.tre')
        self.treeDerepOut = os.path.join(outputDir, 'genome_tree.derep.tre')
        self.treeDerepRootedOut = os.path.join(outputDir, 'genome_tree.derep.rooted.tre')
        self.treeDerepBootstrapOut = os.path.join(outputDir, 'genome_tree.derep.bs.tre')
        self.treeDerepFinalOut = os.path.join(outputDir, 'genome_tree.final.tre')
        self.taxonomyOut = os.path.join(outputDir, 'genome_tree.taxonomy.tsv')
        self.treeMetadata = os.path.join(outputDir, 'genome_tree.metadata.tsv')
        self.phyloHMMsOut = os.path.join(outputDir, 'phylo.hmm')
        self.derepSeqFile = os.path.join(outputDir, 'genome_tree.derep.txt')

        self.phyloUbiquity = 0.90
        self.phyloSingleCopy = 0.90
        self.paralogAcceptPer = 0.01
        self.consistencyAcceptPer = 0.86
        self.consistencyMinTaxa = 20

        # create output directories
        if False:
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
        if False:
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

            sys.exit()

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
        inferGenomeTree.run(self.finalGeneTreeDir, self.alignmentDir, '.aln.masked.faa', self.concatenatedAlignFile, self.treeOut, self.taxonomyOut)

        if False:
            # root genome tree between archaea and bacteria
            print ''
            print '--- Rooting full genome tree ---'
            rerootTree = RerootTree()
            rerootTree.run(self.treeOut, self.treeRootedOut)

            # decorate genome tree with taxonomy using nlevel from tax2tree
            print ''
            print '--- Decorating full genome tree with taxonomic information using tax2tree ---'
            os.system('nlevel -t %s -m %s -o %s' % (self.treeRootedOut, self.taxonomyOut, self.treeTaxonomyOut))

            # dereplicate identical sequences
            print ''
            print '--- Identifying duplicate sequences ---'
            os.system('seqmagick convert --deduplicate-sequences --deduplicated-sequences-file ' + self.derepSeqFile + ' ' + self.concatenatedAlignFile + ' ' + self.derepConcatenatedAlignFile)

            # infer dereplicated genome tree
            print ''
            print '--- Inferring dereplicated genome tree ---'
            outputLog = self.treeDerepOut[0:self.treeDerepOut.rfind('.')] + '.log'
            #cmd = 'FastTreeMP -nosupport -wag -gamma -log ' + outputLog + ' ' + self.derepConcatenatedAlignFile + ' > ' + self.treeDerepOut
            cmd = 'FastTreeMP -wag -gamma -log ' + outputLog + ' ' + self.derepConcatenatedAlignFile + ' > ' + self.treeDerepOut
            os.system(cmd)

            # root genome tree between archaea and bacteria
            print ''
            print '--- Rooting dereplicated genome tree ---'
            rerootTree = RerootTree()
            rerootTree.run(self.treeDerepOut, self.treeDerepRootedOut)

            # calculate bootstraps for genome tree
            print ''
            print '--- Calculating bootstrap support ---'
            #bootstrapTree = BootstrapTree()
            #bootstrapTree.run(self.bootstrapDir, self.treeDerepRootedOut, self.concatenatedAlignFile, 100, numThreads, self.treeDerepBootstrapOut)

            #os.system('cp ' + self.treeDerepBootstrapOut + ' ' + self.treeDerepFinalOut)

            # just use FastTree support values
            os.system('cp ' + self.treeDerepRootedOut + ' ' + self.treeDerepFinalOut)
    
            # decorate dereplicated tree with unique IDs and a complementary file indicating properties of each internal node
            print ''
            print '--- Decorating final tree with lineage-specific statistics and marker set information ---'
            decorateTree = DecorateTree()
            decorateTree.decorate(self.treeTaxonomyOut, self.derepSeqFile, self.treeDerepFinalOut, self.treeMetadata, numThreads)

if __name__ == '__main__':
    print 'GenomeTreeWorkflow v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)

    args = parser.parse_args()

    genomeTreeWorkflow = GenomeTreeWorkflow(args.output_dir)
    genomeTreeWorkflow.run(args.threads)
