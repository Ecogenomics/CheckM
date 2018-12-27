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
import random

import dendropy
from  dendropy.dataobject.taxon import Taxon

from numpy import mean, std, abs, percentile, array

from checkm.lib.img import IMG
from checkm.lib.seqUtils import readFastaBases
from lib.markerSetBuilder import MarkerSetBuilder

class Simulation(object):
    def __init__(self):
        self.markerSetBuilder = MarkerSetBuilder()
        self.img = IMG()
        
        self.simContigLen = 10000
        
    def __selectMarkerSet(self, tree, internalNode, metadata, ubiquityThreshold, singleCopyThreshold, queueOut):
        """Select marker set for parent edge of specified internal node."""
        
        # get genomes descendant from each child of the specified internal node
        leaves = []
        for child in internalNode.child_nodes(): 
            genomeIds = set()  
            for leaf in child.leaf_nodes():
                genomeId = leaf.taxon.label.replace('IMG_', '')
                genomeIds.add(genomeId)
                
                duplicateGenomes = self.markerSetBuilder.duplicateSeqs.get(leaf.taxon.label, [])
                for dup in duplicateGenomes:
                    dupId = dup.replace('IMG_', '')
                    genomeIds.add(dupId)
                 
            leaves.append(genomeIds)
            
        # make sure each set of leaves contains at least a minimum number of genomes
        orderedLeaves = sorted(leaves, key=len)
        if len(orderedLeaves[0]) < 5:
            queueOut.put(('NA', -1, -1, -1, -1, -1))
            return
                   
        # calculate marker genes with all genomes in lineage with the fewest genomes removed 
        binMarkerGenes, _ = self.markerSetBuilder.buildBinMarkerSet(tree, internalNode, ubiquityThreshold, singleCopyThreshold, bMarkerSet = False, genomeIdsToRemove = orderedLeaves[0])
        
        # evaluate accuracy of completeness and contamination estimations on different partial genomes from lineage with fewest genomes   
        testGenomeIds = random.sample(orderedLeaves[0], min(len(orderedLeaves[0]), 100))    
        
        deltaComp = defaultdict(list)
        deltaCont = defaultdict(list)
        
        for testGenomeId in testGenomeIds:   
            geneDistTable = self.img.geneDistTable([testGenomeId], binMarkerGenes.getMarkerGenes(), spacingBetweenContigs=0)
            genomeSize = readFastaBases(os.path.join(self.img.genomeDir, testGenomeId, testGenomeId + '.fna'))
            
            repsPerGenome = 100
            for _ in range(0, repsPerGenome): 
                testComp = random.uniform(0.5, 1.0)
                testCont = random.uniform(0, 0.2)
                trueComp, trueCont, startPartialGenomeContigs = self.markerSetBuilder.sampleGenome(genomeSize, testComp, testCont, self.simContigLen)   
      
                for ms in binMarkerGenes.markerSetIter():  
                    containedMarkerGenes = self.markerSetBuilder.containedMarkerGenes(ms.getMarkerGenes(), geneDistTable[testGenomeId], startPartialGenomeContigs, self.simContigLen)
                    completeness, contamination = ms.genomeCheck(containedMarkerGenes, bIndividualMarkers=True)      
                    if completeness == 0.0:
                        print(ms.getMarkerGenes())
                        print(geneDistTable[testGenomeId])
                        print(startPartialGenomeContigs)
                        print(genomeSize)
                        print('*****************' + testGenomeId)
                        sys.exit()
                    deltaComp[ms.lineageStr].append(completeness - trueComp)
                    deltaCont[ms.lineageStr].append(contamination - trueCont)
            
        # determine lineage-specific marker set with best average performance
        curBest = 1000
        bestUID = None
        dCompBest = 0
        dContBest = 0
        
        for lineageStr in deltaComp:
            dComp, dCont = mean(abs(array(deltaComp[lineageStr]))), mean(abs(array(deltaCont[lineageStr])))

            if (dComp + dCont) < curBest:
                dCompBest = dComp
                dContBest = dCont
                dCompStdBest = std(abs(array(deltaComp[lineageStr])))
                dContStdBest = std(abs(array(deltaCont[lineageStr])))
                bestUID = lineageStr.split('|')[0]
                curBest = dComp + dCont

        queueOut.put((internalNode, bestUID, dCompBest, dCompStdBest, dContBest, dContStdBest))
                        
    def __workerThread(self, tree, metadata, ubiquityThreshold, singleCopyThreshold, queueIn, queueOut):
        """Process each data item in parallel."""

        while True:
            internalNode = queueIn.get(block=True, timeout=None)
            if internalNode == None:
                break
            
            self.__selectMarkerSet(tree, internalNode, metadata, ubiquityThreshold, singleCopyThreshold, queueOut)      
                      
    def __writerThread(self, numInternalNodes, writerQueue):
        """Store or write results of worker threads in a single thread."""

        fout = open('/tmp/simInferBestMarkerSet.tsv', 'w')
        fout.write('Internal node ID\tMarker set ID\tmean % delta comp\tstd % delta comp\tmean % delta cont\tstd % delta cont\n')

        itemsProcessed = 0
        while True:
            internalNode, bestUID, dCompBest, dCompStdBest, dContBest, dContStdBest = writerQueue.get(block=True, timeout=None)
            if internalNode == None:
                break
            
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) internal branches.' % (itemsProcessed, numInternalNodes, float(itemsProcessed)*100/(numInternalNodes))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            if internalNode != 'NA':
                fout.write(internalNode.label + '\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n' % (bestUID, dCompBest, dCompStdBest, dContBest, dContStdBest)) 
            
        fout.close()

        sys.stdout.write('\n')

    def run(self, ubiquityThreshold, singleCopyThreshold, numThreads):
        random.seed(0)
          
        print('\n  Calculating global gene count table.')
        metadata = self.img.genomeMetadata()
        self.markerSetBuilder.globalGeneCountTable = self.img.geneCountTable(list(metadata.keys()))
          
        print('\n  Reading reference genome tree.')
        treeFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
            
        print('  Evaluating %d internal nodes.' % len(tree.internal_nodes()))
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for internalNode in tree.internal_nodes():
            if internalNode.parent_node != None:
                workerQueue.put(internalNode)

        for _ in range(numThreads):
            workerQueue.put(None)

        metadata = self.img.genomeMetadata()
        workerProc = [mp.Process(target = self.__workerThread, args = (tree, metadata, ubiquityThreshold, singleCopyThreshold, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target = self.__writerThread, args = (len(tree.internal_nodes())-1, writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put((None, None, None, None, None, None))
        writeProc.join()
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-t', '--threads', help='Threads to use', type=int, default = 40)

    args = parser.parse_args()

    simulation = Simulation()
    simulation.run(args.ubiquity, args.single_copy, args.threads)
