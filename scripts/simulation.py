#!/usr/bin/env python3

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

"""
Simulate performance of marker sets under different conditions.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict
import gzip
import time
import random

import dendropy
from  dendropy.dataobject.taxon import Taxon

from numpy import mean, std, abs

from checkm.util.img import IMG
from checkm.util.seqUtils import readFastaBases
from lib.markerSetBuilder import MarkerSetBuilder

class Simulation(object):
    def __init__(self):
        self.markerSetBuilder = MarkerSetBuilder()
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')

        self.contigLens = [1000, 2000, 5000, 10000, 20000, 50000]
        self.percentComps = [0.5, 0.7, 0.8, 0.9, 0.95, 1.0]
        self.percentConts = [0.0, 0.05, 0.1, 0.15, 0.2]

    def __workerThread(self, tree, metadata, ubiquityThreshold, singleCopyThreshold, numReplicates, queueIn, queueOut):
        """Process each data item in parallel."""

        while True:
            testGenomeId = queueIn.get(block=True, timeout=None)
            if testGenomeId == None:
                break

            # build marker sets for evaluating test genome
            testNode = tree.find_node_with_taxon_label('IMG_' + testGenomeId)
            binMarkerSets, refinedBinMarkerSet = self.markerSetBuilder.buildBinMarkerSet(tree, testNode.parent_node, ubiquityThreshold, singleCopyThreshold, bMarkerSet=True, genomeIdsToRemove=[testGenomeId])
            #!!!binMarkerSets, refinedBinMarkerSet = self.markerSetBuilder.buildDomainMarkerSet(tree, testNode.parent_node, ubiquityThreshold, singleCopyThreshold, bMarkerSet = False, genomeIdsToRemove = [testGenomeId])

            # determine distribution of all marker genes within the test genome
            geneDistTable = self.img.geneDistTable([testGenomeId], binMarkerSets.getMarkerGenes(), spacingBetweenContigs=0)

            print(('# marker genes: ', len(binMarkerSets.getMarkerGenes())))
            print(('# genes in table: ', len(geneDistTable[testGenomeId])))

            # estimate completeness of unmodified genome
            unmodifiedComp = {}
            unmodifiedCont = {}
            for ms in binMarkerSets.markerSetIter():
                hits = {}
                for mg in ms.getMarkerGenes():
                    if mg in geneDistTable[testGenomeId]:
                        hits[mg] = geneDistTable[testGenomeId][mg]
                completeness, contamination = ms.genomeCheck(hits, bIndividualMarkers=True)
                unmodifiedComp[ms.lineageStr] = completeness
                unmodifiedCont[ms.lineageStr] = contamination

            print((completeness, contamination))

            # estimate completion and contamination of genome after subsampling using both the domain and lineage-specific marker sets
            genomeSize = readFastaBases(os.path.join(self.img.genomeDir, testGenomeId, testGenomeId + '.fna'))
            print(('genomeSize', genomeSize))

            for contigLen in self.contigLens:
                for percentComp in self.percentComps:
                    for percentCont in self.percentConts:
                        deltaComp = defaultdict(list)
                        deltaCont = defaultdict(list)
                        deltaCompSet = defaultdict(list)
                        deltaContSet = defaultdict(list)

                        deltaCompRefined = defaultdict(list)
                        deltaContRefined = defaultdict(list)
                        deltaCompSetRefined = defaultdict(list)
                        deltaContSetRefined = defaultdict(list)

                        trueComps = []
                        trueConts = []

                        numDescendants = {}

                        for _ in range(0, numReplicates):
                            trueComp, trueCont, startPartialGenomeContigs = self.markerSetBuilder.sampleGenome(genomeSize, percentComp, percentCont, contigLen)
                            print((contigLen, trueComp, trueCont, len(startPartialGenomeContigs)))

                            trueComps.append(trueComp)
                            trueConts.append(trueCont)

                            for ms in binMarkerSets.markerSetIter():
                                numDescendants[ms.lineageStr] = ms.numGenomes

                                containedMarkerGenes = self.markerSetBuilder.containedMarkerGenes(ms.getMarkerGenes(), geneDistTable[testGenomeId], startPartialGenomeContigs, contigLen)
                                completeness, contamination = ms.genomeCheck(containedMarkerGenes, bIndividualMarkers=True)
                                deltaComp[ms.lineageStr].append(completeness - trueComp)
                                deltaCont[ms.lineageStr].append(contamination - trueCont)

                                completeness, contamination = ms.genomeCheck(containedMarkerGenes, bIndividualMarkers=False)
                                deltaCompSet[ms.lineageStr].append(completeness - trueComp)
                                deltaContSet[ms.lineageStr].append(contamination - trueCont)

                            for ms in refinedBinMarkerSet.markerSetIter():
                                containedMarkerGenes = self.markerSetBuilder.containedMarkerGenes(ms.getMarkerGenes(), geneDistTable[testGenomeId], startPartialGenomeContigs, contigLen)
                                completeness, contamination = ms.genomeCheck(containedMarkerGenes, bIndividualMarkers=True)
                                deltaCompRefined[ms.lineageStr].append(completeness - trueComp)
                                deltaContRefined[ms.lineageStr].append(contamination - trueCont)

                                completeness, contamination = ms.genomeCheck(containedMarkerGenes, bIndividualMarkers=False)
                                deltaCompSetRefined[ms.lineageStr].append(completeness - trueComp)
                                deltaContSetRefined[ms.lineageStr].append(contamination - trueCont)

                        taxonomy = ';'.join(metadata[testGenomeId]['taxonomy'])
                        queueOut.put((testGenomeId, contigLen, percentComp, percentCont, taxonomy, numDescendants, unmodifiedComp, unmodifiedCont, trueComps, trueConts, deltaComp, deltaCont, deltaCompSet, deltaContSet, deltaCompRefined, deltaContRefined, deltaCompSetRefined, deltaContSetRefined, trueComps, trueConts))

    def __writerThread(self, numTestGenomes, writerQueue):
        """Store or write results of worker threads in a single thread."""

        # summaryOut = open('/tmp/simulation.draft.summary.w_refinement_50.tsv', 'w')
        summaryOut = open('/tmp/simulation.summary.testing.tsv', 'w')
        summaryOut.write('Genome Id\tContig len\t% comp\t% cont')
        summaryOut.write('\tTaxonomy\tMarker set\t# descendants')
        summaryOut.write('\tUnmodified comp\tUnmodified cont\tTrue comp\tTrue cont')
        summaryOut.write('\tIM comp\tIM comp std\tIM cont\tIM cont std')
        summaryOut.write('\tMS comp\tMS comp std\tMS cont\tMS cont std')
        summaryOut.write('\tRIM comp\tRIM comp std\tRIM cont\tRIM cont std')
        summaryOut.write('\tRMS comp\tRMS comp std\tRMS cont\tRMS cont std\n')

        # fout = gzip.open('/tmp/simulation.draft.w_refinement_50.tsv.gz', 'wb')
        fout = gzip.open('/tmp/simulation.testing.tsv.gz', 'wb')
        fout.write('Genome Id\tContig len\t% comp\t% cont')
        fout.write('\tTaxonomy\tMarker set\t# descendants')
        fout.write('\tUnmodified comp\tUnmodified cont\tTrue comp\tTrue cont')
        fout.write('\tIM comp\tIM cont')
        fout.write('\tMS comp\tMS cont')
        fout.write('\tRIM comp\tRIM cont')
        fout.write('\tRMS comp\tRMS cont\tTrue Comp\tTrue Cont\n')

        testsPerGenome = len(self.contigLens) * len(self.percentComps) * len(self.percentConts)

        itemsProcessed = 0
        while True:
            testGenomeId, contigLen, percentComp, percentCont, taxonomy, numDescendants, unmodifiedComp, unmodifiedCont, trueComps, trueConts, deltaComp, deltaCont, deltaCompSet, deltaContSet, deltaCompRefined, deltaContRefined, deltaCompSetRefined, deltaContSetRefined, trueComps, trueConts = writerQueue.get(block=True, timeout=None)
            if testGenomeId == None:
                break

            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, numTestGenomes * testsPerGenome, float(itemsProcessed) * 100 / (numTestGenomes * testsPerGenome))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            for markerSetId in unmodifiedComp:
                summaryOut.write(testGenomeId + '\t%d\t%.2f\t%.2f' % (contigLen, percentComp, percentCont))
                summaryOut.write('\t' + taxonomy + '\t' + markerSetId + '\t' + str(numDescendants[markerSetId]))
                summaryOut.write('\t%.3f\t%.3f' % (unmodifiedComp[markerSetId], unmodifiedCont[markerSetId]))
                summaryOut.write('\t%.3f\t%.3f' % (mean(trueComps), std(trueConts)))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaComp[markerSetId])), std(abs(deltaComp[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaCont[markerSetId])), std(abs(deltaCont[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaCompSet[markerSetId])), std(abs(deltaCompSet[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaContSet[markerSetId])), std(abs(deltaContSet[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaCompRefined[markerSetId])), std(abs(deltaCompRefined[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaContRefined[markerSetId])), std(abs(deltaContRefined[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaCompSetRefined[markerSetId])), std(abs(deltaCompSetRefined[markerSetId]))))
                summaryOut.write('\t%.3f\t%.3f' % (mean(abs(deltaContSetRefined[markerSetId])), std(abs(deltaContSetRefined[markerSetId]))))
                summaryOut.write('\n')

                fout.write(testGenomeId + '\t%d\t%.2f\t%.2f' % (contigLen, percentComp, percentCont))
                fout.write('\t' + taxonomy + '\t' + markerSetId + '\t' + str(numDescendants[markerSetId]))
                fout.write('\t%.3f\t%.3f' % (unmodifiedComp[markerSetId], unmodifiedCont[markerSetId]))
                fout.write('\t%s' % ','.join(map(str, trueComps)))
                fout.write('\t%s' % ','.join(map(str, trueConts)))
                fout.write('\t%s' % ','.join(map(str, deltaComp[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaCont[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaCompSet[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaContSet[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaCompRefined[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaContRefined[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaCompSetRefined[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, deltaContSetRefined[markerSetId])))
                fout.write('\t%s' % ','.join(map(str, trueComps)))
                fout.write('\t%s' % ','.join(map(str, trueConts)))
                fout.write('\n')

        summaryOut.close()
        fout.close()

        sys.stdout.write('\n')

    def run(self, ubiquityThreshold, singleCopyThreshold, numReplicates, numThreads):
        print('\n  Reading reference genome tree.')
        treeFile = os.path.join('/srv', 'db', 'checkm', 'genome_tree', 'genome_tree_full.refpkg', 'genome_tree.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)

        print(('    Number of taxa in tree: %d' % (len(tree.leaf_nodes()))))

        genomesInTree = set()
        for leaf in tree.leaf_iter():
            genomesInTree.add(leaf.taxon.label.replace('IMG_', ''))

        # get all draft genomes for testing
        print('')
        metadata = self.img.genomeMetadata()
        print(('  Total genomes: %d' % len(metadata)))

        genomeIdsToTest = genomesInTree - self.img.filterGenomeIds(genomesInTree, metadata, 'status', 'Finished')
        print(('  Number of draft genomes: %d' % len(genomeIdsToTest)))

        print('')
        print('  Pre-computing genome information for calculating marker sets:')
        start = time.time()
        self.markerSetBuilder.readLineageSpecificGenesToRemove()
        end = time.time()
        print(('    readLineageSpecificGenesToRemove: %.2f' % (end - start)))


        start = time.time()
        # self.markerSetBuilder.cachedGeneCountTable = self.img.geneCountTable(metadata.keys())
        end = time.time()
        print(('    globalGeneCountTable: %.2f' % (end - start)))

        start = time.time()
        # self.markerSetBuilder.precomputeGenomeSeqLens(metadata.keys())
        end = time.time()
        print(('    precomputeGenomeSeqLens: %.2f' % (end - start)))

        start = time.time()
        # self.markerSetBuilder.precomputeGenomeFamilyPositions(metadata.keys(), 0)
        end = time.time()
        print(('    precomputeGenomeFamilyPositions: %.2f' % (end - start)))

        print('')
        print(('  Evaluating %d test genomes.' % len(genomeIdsToTest)))
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for testGenomeId in genomeIdsToTest:
            workerQueue.put(testGenomeId)

        for _ in range(numThreads):
            workerQueue.put(None)

        workerProc = [mp.Process(target=self.__workerThread, args=(tree, metadata, ubiquityThreshold, singleCopyThreshold, numReplicates, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target=self.__writerThread, args=(len(genomeIdsToTest), writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put((None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None))
        writeProc.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default=0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default=0.97)
    parser.add_argument('-x', '--replicates', help='Replicates per genome.', type=int, default=20)
    parser.add_argument('-t', '--threads', help='Threads to use', type=int, default=46)

    args = parser.parse_args()

    simulation = Simulation()
    simulation.run(args.ubiquity, args.single_copy, args.replicates, args.threads)
