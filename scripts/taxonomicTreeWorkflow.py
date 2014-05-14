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

__prog_desc__ = 'complete workflow for inferring a tree for a taxonomic group'

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
import multiprocessing as mp
import random

from checkm.lib.img import IMG

from lib.markerSetBuilder import MarkerSetBuilder

from checkm.seqUtils import readFasta
from checkm.hmmer import HMMERRunner

from genometreeworkflow.makeTrees import MakeTrees
from genometreeworkflow.paralogTest import ParalogTest
from genometreeworkflow.consistencyTest import ConsistencyTest
from genometreeworkflow.getPhylogeneticHMMs import GetPhylogeneticHMMs
from genometreeworkflow.inferGenomeTree import InferGenomeTree

class GenomeTreeWorkflow(object):
    def __init__(self, outputDir):
        self.img = IMG()
        self.markerSetBuilder = MarkerSetBuilder()

        if os.path.exists(outputDir):
            print '[Error] Output directory already exists: ' + outputDir
            sys.exit(0)
        else:
            os.makedirs(outputDir)

        self.__checkForHMMER()
        self.__checkForFastTree()

        self.hmmDir = os.path.join(outputDir, 'phylo_hmms')
        self.alignmentDir = os.path.join(outputDir, 'gene_alignments')
        self.geneTreeDir = os.path.join(outputDir, 'gene_trees')
        self.conspecificGeneTreeDir = os.path.join(outputDir, 'gene_trees_conspecific')
        self.finalGeneTreeDir = os.path.join(outputDir, 'gene_trees_final')

        self.consistencyOut = os.path.join(outputDir, 'genome_tree.consistency.tsv')
        self.concatenatedAlignFile = os.path.join(outputDir, 'genome_tree.concatenated.faa')
        self.derepConcatenatedAlignFile = os.path.join(outputDir, 'genome_tree.concatenated.derep.fasta')
        self.treeOut = os.path.join(outputDir, 'genome_tree.tre')
        self.treeOutAce = os.path.join(outputDir, 'genome_tree.ace_ids.tre')

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
        #self.consistencyAcceptPer = 0.95    # for trees at the class-level
        self.consistencyAcceptPer = 0.906   # for trees at the phylum-level
        self.consistencyMinTaxa = 20

        # create output directories
        os.makedirs(self.hmmDir)
        os.makedirs(self.alignmentDir)
        os.makedirs(self.geneTreeDir)
        os.makedirs(self.conspecificGeneTreeDir)
        os.makedirs(self.finalGeneTreeDir)

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

    def __genesInGenomes(self, genomeIds):
        genesInGenomes = {}
        for genomeId in genomeIds:
            markerIdToGeneIds = defaultdict(set)
            for line in open(os.path.join(IMG.genomeDir, genomeId, genomeId + IMG.pfamExtension)):
                lineSplit = line.split('\t')
                markerIdToGeneIds[lineSplit[8]].add(lineSplit[0])

            for line in open(os.path.join(IMG.genomeDir, genomeId, genomeId + IMG.tigrExtension)):
                lineSplit = line.split('\t')
                markerIdToGeneIds[lineSplit[6]].add(lineSplit[0])

            genesInGenomes[genomeId] = markerIdToGeneIds

        return genesInGenomes

    def __fetchMarkerModels(self, universalMarkerGenes, outputModelDir):
        markerIdToName = {}
        for line in open('/srv/whitlam/bio/db/pfam/27/Pfam-A.hmm'):
            if 'NAME' in line:
                name = line.split()[1].rstrip()
            elif 'ACC' in line:
                acc = line.split()[1].rstrip()
                markerId = acc.replace('PF', 'pfam')
                markerId = markerId[0:markerId.rfind('.')]
                markerIdToName[markerId] = name

        for markerId in universalMarkerGenes:
            if 'pfam' in markerId:
                os.system('hmmfetch  /srv/whitlam/bio/db/pfam/27/Pfam-A.hmm ' + markerIdToName[markerId] + ' > ' + os.path.join(outputModelDir, markerId.replace('pfam', 'PF') + '.hmm'))
            else:
                os.system('hmmfetch  /srv/whitlam/bio/db/tigrfam/13.0/' + markerId + '.HMM ' + markerId + ' > ' + os.path.join(outputModelDir, markerId + '.hmm'))

    def __alignMarkers(self, genomeIds, markerGenes, genesInGenomes, numThreads, outputGeneDir, outputModelDir):
        """Perform multithreaded alignment of marker genes using HMM align."""

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for _, markerId in enumerate(markerGenes):
            workerQueue.put(markerId)

        for _ in range(numThreads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__runHmmAlign, args = (genomeIds, genesInGenomes, outputGeneDir, outputModelDir, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target = self.__reportThreads, args = (len(markerGenes), writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put(None)
        writeProc.join()

    def __runHmmAlign(self, genomeIds, genesInGenomes, outputGeneDir, outputModelDir, queueIn, queueOut):
        """Run each marker gene in a separate thread."""

        while True:
            markerId = queueIn.get(block=True, timeout=None)
            if markerId == None:
                break

            modelName = markerId
            if modelName.startswith('pfam'):
                modelName = modelName.replace('pfam', 'PF')

            markerSeqFile = os.path.join(outputGeneDir, modelName + '.faa')
            fout = open(markerSeqFile, 'w')
            for genomeId in genomeIds:
                seqs = readFasta(IMG.genomeDir + '/' + genomeId + '/' + genomeId + '.genes.faa')

                for geneId in genesInGenomes[genomeId].get(markerId, []):
                    if geneId not in seqs:
                        # this shouldn't be necessary, but the IMG metadata isn't always
                        # perfectly in sync with the sequence data
                        continue

                    fout.write('>' + genomeId + '|' + geneId + '\n')
                    fout.write(seqs[geneId] + '\n')
            fout.close()

            hmmer = HMMERRunner('align')
            hmmer.align(os.path.join(outputModelDir, modelName + '.hmm'), markerSeqFile, os.path.join(outputGeneDir, modelName + '.aln.faa'), trim=False, outputFormat='Pfam')
            self.__maskAlignment(os.path.join(outputGeneDir, modelName + '.aln.faa'), os.path.join(outputGeneDir, modelName + '.aln.masked.faa'))

            queueOut.put(modelName)

    def __reportThreads(self, numGenes, writerQueue):
        """Store confidence intervals (i.e., to shared memory)."""

        numProcessedGenes = 0
        while True:
            markerId = writerQueue.get(block=True, timeout=None)
            if markerId == None:
                break

            numProcessedGenes += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) marker genes.' % (numProcessedGenes, numGenes, float(numProcessedGenes)*100/numGenes)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def __maskAlignment(self, inputFile, outputFile):
        """Read HMMER alignment in STOCKHOLM format and output masked alignment in FASTA format."""
        # read STOCKHOLM alignment
        seqs = {}
        for line in open(inputFile):
            line = line.rstrip()
            if line == '' or line[0] == '#' or line == '//':
                if 'GC RF' in line:
                    mask = line.split('GC RF')[1].strip()
                continue
            else:
                lineSplit = line.split()
                seqs[lineSplit[0]] = lineSplit[1].upper().replace('.', '-').strip()

        # output masked sequences in FASTA format
        fout = open(outputFile, 'w')
        for seqId, seq in seqs.iteritems():
            fout.write('>' + seqId + '\n')

            maskedSeq = ''.join([seq[i] for i in xrange(0, len(seq)) if mask[i] == 'x'])
            fout.write(maskedSeq + '\n')
        fout.close()

    def __taxonomicMarkers(self, taxaStr, metadata, ubiquityThreshold, singleCopyThreshold):
        """Get genomes and marker genes for a specific lineage."""

        genomeIds = self.img.genomeIdsByTaxonomy(taxaStr, metadata)
        
        # hack to add in other genomes with incorrect taxonomy
        aceIds = set()
        for line in open('./taxonomicTrees/firmicutes.nds'):
            aceIds.add(line.strip())
            
        aceIdsToImgIds = self.aceIdsToImgIds()
        for aceId in aceIds:
            if aceId in aceIdsToImgIds:
                genomeId = aceIdsToImgIds[aceId]
                if genomeId in metadata:
                    genomeIds.add(genomeId)
        
        markerGenes = self.markerSetBuilder.buildMarkerGenes(genomeIds, ubiquityThreshold, singleCopyThreshold )

        return genomeIds, markerGenes

    def inferGeneTrees(self, phyloUbiquityThreshold, phyloSingleCopyThreshold, numThreads, outputGeneDir, outputModelDir, outgroupSize):
        # make sure output directory is empty
        if not os.path.exists(outputGeneDir):
            os.makedirs(outputGeneDir)

        if not os.path.exists(outputModelDir):
            os.makedirs(outputModelDir)

        files = os.listdir(outputGeneDir)
        for f in files:
            os.remove(os.path.join(outputGeneDir, f))

        # get genomes and marker genes for taxonomic groups of interest
        print ''
        print 'Identifying genomes and marker genes of interest:'
        metadata = self.img.genomeMetadata()
        ingroupGenomeIds, ingroupMarkers = self.__taxonomicMarkers('Bacteria;Firmicutes', metadata, phyloUbiquityThreshold, phyloSingleCopyThreshold)
        outgroupGenomeIds = self.img.genomeIdsByTaxonomy('Bacteria;Coprothermobacter', metadata)
        #alphaGenomeIds, _ = self.__taxonomicMarkers('Bacteria;Proteobacteria;Alphaproteobacteria', metadata, phyloUbiquityThreshold, phyloSingleCopyThreshold)

        print '  Identified ingroup genomes: %d' % len(ingroupGenomeIds)
        print '  Identified outgroup genomes: %d' % len(outgroupGenomeIds)

        numOutgroupTaxa = min(outgroupSize, len(outgroupGenomeIds))
        print ''
        print '  Selecting %d taxa from the outgroup.' % (numOutgroupTaxa)
        genomeIds = ingroupGenomeIds.union(random.sample(outgroupGenomeIds, numOutgroupTaxa))

        self.imgIdsToAceIds(genomeIds)

        print '  Identified markers: %d' % len(ingroupMarkers)

        # get mapping of marker ids to gene ids for each genome
        print '  Determine genes for genomes of interest.'
        genesInGenomes = self.__genesInGenomes(genomeIds)

        # get HMM for each marker gene
        print '  Fetching HMM for each marker genes.'
        self.__fetchMarkerModels(ingroupMarkers, outputModelDir)

        # align gene sequences and infer gene trees
        print '  Aligning marker genes:'
        self.__alignMarkers(genomeIds, ingroupMarkers, genesInGenomes, numThreads, outputGeneDir, outputModelDir)

        return genomeIds

    def imgIdsToAceIds(self, imgIds):
        imgIdToAceId = {}
        for line in open('ggg_tax_img.feb_2014.txt'):
            lineSplit = line.split('\t')
            imgIdToAceId[lineSplit[1].rstrip()] = lineSplit[0]

        missing = 0
        for imgId in imgIds:
            if imgId not in imgIdToAceId:
                missing += 1

        print '  Number of genomes without an ACE id: ' + str(missing)

        return imgIdToAceId
    
    def aceIdsToImgIds(self):
        aceIdsToImgIds = {}
        for line in open('ggg_tax_img.feb_2014.txt'):
            lineSplit = line.split('\t')
            aceIdsToImgIds[lineSplit[0].strip()] = lineSplit[1].strip()

        return aceIdsToImgIds

    def run(self, numThreads, outgroupSize):

        # identify genes suitable for phylogenetic inference
        print '--- Identifying genes suitable for phylogenetic inference ---'
        genomeIds = self.inferGeneTrees(self.phyloUbiquity, self.phyloSingleCopy, numThreads, self.alignmentDir, self.hmmDir, outgroupSize)

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
        inferGenomeTree.run(self.finalGeneTreeDir, self.alignmentDir, '.aln.masked.faa', self.concatenatedAlignFile, self.treeOut, self.taxonomyOut, bSupportValues = True)

        # replace IMG identifiers with ACE identifiers
        imgIdToAceId = self.imgIdsToAceIds(genomeIds)
        with open(self.treeOut) as f:
            tree = ''.join(f.readlines())

            for genomeId in genomeIds:
                if genomeId in imgIdToAceId:
                    tree = tree.replace('IMG_' + genomeId, imgIdToAceId[genomeId])

        fout = open(self.treeOutAce, 'w')
        fout.write(tree)
        fout.close()

if __name__ == '__main__':
    print 'GenomeTreeWorkflow v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)
    parser.add_argument('-n', '--outgroup_size', help='number of taxa in outgroup', type = int, default = 30)

    args = parser.parse_args()

    genomeTreeWorkflow = GenomeTreeWorkflow(args.output_dir)
    genomeTreeWorkflow.run(args.threads, args.outgroup_size)
