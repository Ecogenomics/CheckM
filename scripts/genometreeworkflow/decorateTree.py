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

__prog_desc__ = 'decorating tree with lineage-specific statistics and marker set information'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0. 0.2'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse
import multiprocessing as mp

import dendropy

from numpy import mean, std

from checkm.util.img import IMG
from .markerSetBuilder import MarkerSetBuilder

class DecorateTree(object):
    def __init__(self):
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        self.pfamHMMs = '/srv/whitlam/bio/db/pfam/27/Pfam-A.hmm'
        self.markerSetBuilder = MarkerSetBuilder()

    def __meanStd(self, metadata, genomeIds, category):
        values = []
        for genomeId in genomeIds:
            genomeId = genomeId.replace('IMG_', '')
            v = metadata[genomeId][category]
            if v != 'NA':
                values.append(v)

        return mean(values), std(values)

    def __calculateMarkerSet(self, genomeLabels, ubiquityThreshold=0.97, singleCopyThreshold=0.97):
        """Calculate marker set for a set of genomes."""

        # get genome IDs from genome labels
        genomeIds = set()
        for genomeLabel in genomeLabels:
            genomeIds.add(genomeLabel.replace('IMG_', ''))

        markerSet = self.markerSetBuilder.buildMarkerSet(genomeIds, ubiquityThreshold, singleCopyThreshold)

        return markerSet.markerSet

    def __pfamIdToPfamAcc(self, img):
        pfamIdToPfamAcc = {}
        for line in open(self.pfamHMMs):
            if 'ACC' in line:
                acc = line.split()[1].strip()
                pfamId = acc.split('.')[0]

                pfamIdToPfamAcc[pfamId] = acc

        return pfamIdToPfamAcc

    def decorate(self, taxaTreeFile, derepFile, inputTreeFile, metadataOut, numThreads):
        # read genome metadata
        print('  Reading metadata.')
        metadata = self.img.genomeMetadata()

        # read list of taxa with duplicate sequences
        print('  Read list of taxa with duplicate sequences.')
        duplicateTaxa = {}
        for line in open(derepFile):
            lineSplit = line.rstrip().split()
            if len(lineSplit) > 1:
                duplicateTaxa[lineSplit[0]] = lineSplit[1:]

        # build gene count table
        print('  Building gene count table.')
        genomeIds = list(self.img.genomeMetadata().keys())
        print('    # trusted genomes = ' + str(len(genomeIds)))

        # calculate statistics for each internal node using multiple threads
        print('  Calculating statistics for each internal node.')
        self.__internalNodeStatistics(taxaTreeFile, inputTreeFile, duplicateTaxa, metadata, metadataOut, numThreads)

    def __internalNodeStatistics(self, taxaTreeFile, inputTreeFile, duplicateTaxa, metadata, metadataOut, numThreads):

        # determine HMM model accession numbers
        pfamIdToPfamAcc = self.__pfamIdToPfamAcc(self.img)

        taxaTree = dendropy.Tree.get_from_path(taxaTreeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        inputTree = dendropy.Tree.get_from_path(inputTreeFile, schema='newick', as_rooted=True, preserve_underscores=True)

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        uniqueId = 0
        for node in inputTree.internal_nodes():
            uniqueId += 1
            workerQueue.put((uniqueId, node))

        for _ in range(numThreads):
            workerQueue.put((None, None))

        calcProc = [mp.Process(target=self.__processInternalNode, args=(taxaTree, duplicateTaxa, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target=self.__reportStatistics, args=(metadata, metadataOut, inputTree, inputTreeFile, pfamIdToPfamAcc, writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None, None, None, None, None, None))
        writeProc.join()

    def __processInternalNode(self, taxaTree, duplicateTaxa, queueIn, queueOut):
        """Run each marker gene in a separate thread."""

        while True:
            uniqueId, node = queueIn.get(block=True, timeout=None)
            if uniqueId == None:
                break

            # find corresponding internal node in taxa tree
            labels = []
            for leaf in node.leaf_nodes():
                labels.append(leaf.taxon.label)
                if leaf.taxon.label in duplicateTaxa:
                    for genomeId in duplicateTaxa[leaf.taxon.label]:
                        labels.append(genomeId)

            # check if there is a taxonomic label
            mrca = taxaTree.mrca(taxon_labels=labels)
            taxaStr = ''
            if mrca.label:
                taxaStr = mrca.label.replace(' ', '')

            # give node a unique Id while retraining bootstrap value
            bootstrap = ''
            if node.label:
                bootstrap = node.label
            nodeLabel = 'UID' + str(uniqueId) + '|' + taxaStr + '|' + bootstrap

            # calculate marker set
            markerSet = self.__calculateMarkerSet(labels)

            queueOut.put((uniqueId, labels, markerSet, taxaStr, bootstrap, node.oid, nodeLabel))

    def __reportStatistics(self, metadata, metadataOut, inputTree, inputTreeFile, pfamIdToPfamAcc, writerQueue):
        """Store statistics for internal node."""

        fout = open(metadataOut, 'w')
        fout.write('UID\t# genomes\tTaxonomy\tBootstrap')
        fout.write('\tGC mean\tGC std')
        fout.write('\tGenome size mean\tGenome size std')
        fout.write('\tGene count mean\tGene count std')
        fout.write('\tMarker set')
        fout.write('\n')

        numProcessedNodes = 0
        numInternalNodes = len(inputTree.internal_nodes())
        while True:
            uniqueId, labels, markerSet, taxaStr, bootstrap, nodeID, nodeLabel = writerQueue.get(block=True, timeout=None)
            if uniqueId == None:
                break

            numProcessedNodes += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) internal nodes.' % (numProcessedNodes, numInternalNodes, float(numProcessedNodes) * 100 / numInternalNodes)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            fout.write('UID' + str(uniqueId) + '\t' + str(len(labels)) + '\t' + taxaStr + '\t' + bootstrap)

            m, s = self.__meanStd(metadata, labels, 'GC %')
            fout.write('\t' + str(m * 100) + '\t' + str(s * 100))

            m, s = self.__meanStd(metadata, labels, 'genome size')
            fout.write('\t' + str(m) + '\t' + str(s))

            m, s = self.__meanStd(metadata, labels, 'gene count')
            fout.write('\t' + str(m) + '\t' + str(s))

            # change model names to accession numbers, and make
            # sure there is an HMM model for each PFAM
            mungedMarkerSets = []
            for geneSet in markerSet:
                s = set()
                for geneId in geneSet:
                    if 'pfam' in geneId:
                        pfamId = geneId.replace('pfam', 'PF')
                        if pfamId in pfamIdToPfamAcc:
                            s.add(pfamIdToPfamAcc[pfamId])
                    else:
                        s.add(geneId)
                mungedMarkerSets.append(s)

            fout.write('\t' + str(mungedMarkerSets))

            fout.write('\n')

            node = inputTree.find_node(filter_fn=lambda n: hasattr(n, 'oid') and n.oid == nodeID)
            node.label = nodeLabel

        sys.stdout.write('\n')

        fout.close()

        inputTree.write_to_path(inputTreeFile, schema='newick', suppress_rooting=True, unquoted_underscores=True)

if __name__ == '__main__':
    print('DecorateTree v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('taxa_tree', help='tree with internal nodes labeled with taxonomic information')
    parser.add_argument('derep_file', help='file indicating dereplicated nodes in input tree')
    parser.add_argument('input_tree', help='tree to decorate with new labels; assumed to contain bootstrap values')
    parser.add_argument('metadata_file', help='file to contain metadata for each internal node')
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=16)

    args = parser.parse_args()

    decorateTree = DecorateTree()
    decorateTree.decorate(args.taxa_tree, args.derep_file, args.input_tree, args.metadata_file, args.threads)
