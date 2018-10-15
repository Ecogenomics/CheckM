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
Assess performance of marker set selection criteria.
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

import dendropy
from  dendropy.dataobject.taxon import Taxon

from numpy import mean, std, abs, array, percentile

from checkm.lib.img import IMG
from lib.markerSetBuilder import MarkerSetBuilder

class MarkerSetSelection(object):
    def __init__(self):
        self.simFile = './experiments/simulation.tuning.genus.summary.tsv'
        self.looRank = 5
        
        self.markerSetBuilder = MarkerSetBuilder()
        self.img = IMG()
        
    def __stabilityTest(self, genomeIds, ubiquityThreshold = 0.97, singleCopyThreshold = 0.97, stabilityThreshold = 0.05):
        """Test stability of marker set for a group of genomes using LOO-testing."""
        
        # quick escape for lineage that are clearly stable
        if len(genomeIds) > 200:
            return True
        
        # calculate marker sets using a LOO-testing
        looMarkerGenes = []
        for genomeId in genomeIds:
            looGenomeIds = genomeIds.difference([genomeId])
            
            # calculate marker genes
            geneCountTable = self.img.geneCountTable(looGenomeIds)
            markerGenes = self.markerSetBuilder.markerGenes(looGenomeIds, geneCountTable, ubiquityThreshold*len(looGenomeIds), singleCopyThreshold*len(looGenomeIds))
            tigrToRemove = self.img.identifyRedundantTIGRFAMs(markerGenes)
            markerGenes = markerGenes - tigrToRemove
            
            looMarkerGenes.append(markerGenes)
            
        # calculate change in marker set for all pairs
        markerSetSize = []
        diffMarkerSet = []
        for i in range(0, len(looMarkerGenes)):
            markerSetSize.append(len(looMarkerGenes[i]))
            for j in range(i+1, len(looMarkerGenes)):     
                symmDiff = looMarkerGenes[i].symmetric_difference(looMarkerGenes[j])
                diffMarkerSet.append(len(symmDiff))
                            
        print((len(genomeIds), mean(diffMarkerSet), mean(markerSetSize)))
        return (float(mean(diffMarkerSet)) / mean(markerSetSize)) <= stabilityThreshold
        
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
        
    def __distToNodePercentileTest(self, genomeNode, markerSetNode, leaves, percentileTest):
        
        distToBin = self.__distanceToAncestor(genomeNode, markerSetNode)
        
        distToLeaves = []
        for leaf in leaves:
            distToLeaves.append(self.__distanceToAncestor(leaf, markerSetNode))          
               
        return distToBin < percentile(distToLeaves, percentileTest)
     
    def __selectMarkerSetNode(self, tree, genomeId, metadata, taxonToGenomeIds):
        """Determine lineage-specific marker set to use for assessing the giving genome."""
        
        # read genomes removed from tree as a result of duplicate sequences
        duplicateSeqs = self.markerSetBuilder.readDuplicateSeqs()
        
        # determine location of genome in tree     
        node = tree.find_node_with_taxon_label('IMG_' + genomeId)

        # ascend tree to root looking for suitable marker set
        curNode = node.parent_node
        while curNode != None:
            uniqueId = curNode.label.split('|')[0]
            
            genomeIds = set()
            for leaf in curNode.leaf_nodes():
                genomeIds.add(leaf.taxon.label.replace('IMG_', ''))
                
                duplicateGenomes = duplicateSeqs.get(leaf.taxon.label, [])
                for dup in duplicateGenomes:
                    genomeIds.add(dup.replace('IMG_', ''))
                          
            # remove genome (LOO-style analysis)
            print(('Full:', len(genomeIds)))
            genomeIds.difference_update([genomeId])
            print(('LOO:', len(genomeIds)))
            
            # remove all genomes from the same taxonomic group as the genome of interest
            taxon = metadata[genomeId]['taxonomy'][self.looRank]
            genomeIds.difference_update(taxonToGenomeIds[taxon]) 
            print(('Rank reduced:', len(genomeIds)))
              
            print(uniqueId)
            if len(genomeIds) > 10 and self.__stabilityTest(genomeIds):
                uidSelected = uniqueId
                break
                
            curNode = curNode.parent_node
            if curNode == None:
                # reach root so use universal marker set
                uidSelected = uniqueId
            
        return uidSelected
    
    def __bestMarkerSet(self, genomeId, simResults):
        """Get stats for best marker set."""
        curBest = 1000
        bestUID = None
        for uid, results in list(simResults[genomeId].items()):
            numDescendants, dComp, dCont = results
            if (dComp + dCont) < curBest:
                numDescendantsBest = numDescendants
                dCompBest = dComp
                dContBest = dCont
                bestUID = uid
                curBest = dComp + dCont
                
        return bestUID, numDescendantsBest, dCompBest, dContBest
        
    
    def __workerThread(self, tree, simResults, metadata, taxonToGenomeIds, queueIn, queueOut):
        """Process each data item in parallel."""

        while True:
            testGenomeId = queueIn.get(block=True, timeout=None)
            if testGenomeId == None:
                break
            
            uidSelected = self.__selectMarkerSetNode(tree, testGenomeId, metadata, taxonToGenomeIds)
            numDescendantsSelected, dCompSelected, dContSelected = simResults[testGenomeId][uidSelected]
            
            # find best marker set
            bestUID, numDescendantsBest, dCompBest, dContBest = self.__bestMarkerSet(testGenomeId, simResults)

            queueOut.put((testGenomeId, uidSelected, numDescendantsSelected, dCompSelected, dContSelected, bestUID, numDescendantsBest, dCompBest, dContBest))
                      
    def __writerThread(self, numTestGenomes, writerQueue):
        """Store or write results of worker threads in a single thread."""
        
        fout = open('./experiments/markerSetSelection.tsv', 'w')
        
        fout.write('Genome Id\tSelected UID\t# descendants\tSelected dComp\tSelected dCont\tBest UID\t# descendants\tBest dComp\tBest dCont\tdDescendants\tdComp\tdCont\n')
        
        itemsToProcess = 0
        
        dComps = []
        dConts = []
        
        dCompsPer = []
        dContsPer = []
        
        bestComp = []
        bestCont = []
        
        selectedComp = []
        selectedCont = []
        
        while True:
            testGenomeId, uidSelected, numDescendantsSelected, dCompSelected, dContSelected, bestUID, numDescendantsBest, dCompBest, dContBest = writerQueue.get(block=True, timeout=None)
            if testGenomeId == None:
                break

            itemsToProcess += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test genomes.' % (itemsToProcess, numTestGenomes, float(itemsToProcess)*100/(numTestGenomes))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            dComp = abs(dCompSelected - dCompBest)
            dCont = abs(dContSelected - dContBest)
            dDescendants = abs(numDescendantsSelected - numDescendantsBest)
            fout.write('%s\t%s\t%d\t%.4f\t%.4f\t%s\t%d\t%.4f\t%.4f\t%d\t%.4f\t%.4f\n' % (testGenomeId, uidSelected, numDescendantsSelected, dCompSelected, dContSelected, bestUID, numDescendantsBest, dCompBest, dContBest, dDescendants, dComp, dCont))

            dComps.append(dComp)
            dConts.append(dCont)
            
            dCompsPer.append(dComp*100.0 / dCompBest)
            dContsPer.append(dCont*100.0 / max(dContBest, 0.01))
            
            bestComp.append(dCompBest)
            bestCont.append(dContBest)
            
            selectedComp.append(dCompSelected)
            selectedCont.append(dContSelected)

        sys.stdout.write('\n')
        fout.close()
        
        print('')
        print('  General results:')
        print(('   Best comp: %.2f +/- %.2f' % (mean(bestComp), std(bestComp))))
        print(('   Best cont: %.2f +/- %.2f' % (mean(bestCont), std(bestCont))))
        print(('   Selected comp: %.2f +/- %.2f' % (mean(selectedComp), std(selectedComp))))
        print(('   Selected cont: %.2f +/- %.2f' % (mean(selectedCont), std(selectedCont))))
        print('')
        print(('   Delta comp: %.2f +/- %.2f' % (mean(dComps), std(dComps))))
        print(('   Delta cont: %.2f +/- %.2f' % (mean(dConts), std(dConts))))
        print(('   Delta comp per error: %.1f +/- %.1f' % (mean(dCompsPer), std(dCompsPer))))
        print(('   Delta cont per error: %.1f +/- %.1f' % (mean(dContsPer), std(dContsPer))))
        
    def __distanceToAncestor(self, leaf, ancestor):
        dist = 0
        
        curNode = leaf
        while curNode != ancestor:
            dist += curNode.edge_length
            
            curNode = curNode.parent_node
        
        return dist
        
    def __bestNodeProperties(self, genomeId, tree, bestUID):
        # determine location of genome in tree     
        node = tree.find_node_with_taxon_label('IMG_' + genomeId)

        # find node of best marker set
        curNode = node.parent_node
        nodesToBin = 0
        distanceToBin = node.edge_length
        distanceToLeaves = []
        while curNode != None:
            uniqueId = curNode.label.split('|')[0]
            
            nodesToBin += 1  
            
            if uniqueId == bestUID:
                for leaf in curNode.leaf_nodes():
                    if leaf != node:
                        dist = self.__distanceToAncestor(leaf, curNode)
                        distanceToLeaves.append(dist)
                break
            
            distanceToBin += curNode.edge_length
            
            curNode = curNode.parent_node
            
        return nodesToBin, distanceToBin, mean(distanceToLeaves)
        
    def __propertiesOfBestMarkerSets(self, tree, simResults):
        
        numDescendants = []
        nodesToBin = []
        distanceToBin = []
        avgDistanceToLeaf = []
        percDiffs = []
        for genomeId in simResults:
            bestUID, numDescendantsBest, _, _ = self.__bestMarkerSet(genomeId, simResults)
            nodesToBinBest, distanceToBinBest, avgDistanceToLeafBest = self.__bestNodeProperties(genomeId, tree, bestUID)
            
            numDescendants.append(numDescendantsBest)
            nodesToBin.append(nodesToBinBest)
            distanceToBin.append(distanceToBinBest)
            avgDistanceToLeaf.append(avgDistanceToLeafBest)
            
            percDiff = abs(distanceToBinBest - avgDistanceToLeafBest) * 100 / distanceToBinBest
            percDiffs.append(percDiff)
            
        print(('    # descendants: %.2f +/- %.2f' % (mean(numDescendants), std(numDescendants))))
        print(('    # nodes to bin: %.2f +/- %.2f' % (mean(nodesToBin), std(nodesToBin))))
        print(('    Distance to bin: %.2f +/- %.2f' % (mean(distanceToBin), std(distanceToBin))))
        
        distanceToBin = array(distanceToBin)
        avgDistanceToLeaf = array(avgDistanceToLeaf)
        print(('    Distance to bin - average distance to leaf: %.2f +/- %.2f' % (mean(abs(distanceToBin - avgDistanceToLeaf)), std(abs(distanceToBin - avgDistanceToLeaf)))))
        print(('    Percent difference to average leaf distance: %.2f +/- %.2f' % (mean(percDiffs), std(percDiffs))))
        print('')

    def run(self, numThreads):
        # read reference tree
        print('\n  Reading reference genome tree.')
        treeFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # get all genomes with a given taxon label
        metadata = self.img.genomeMetadata()
        taxonToGenomeIds = defaultdict(set)
        for genomeId in metadata:
            for t in metadata[genomeId]['taxonomy']:
                taxonToGenomeIds[t].add(genomeId)

        # read simulation results
        print('  Reading simulation results.')
        
        simResults = defaultdict(dict)
        with open(self.simFile) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                simId = lineSplit[0] + '-' + lineSplit[1] + '-' + lineSplit[2] + '-' + lineSplit[3]
                uid = lineSplit[5].split('|')[0].strip()
                numDescendants = int(lineSplit[6])
                comp = float(lineSplit[21])
                cont = float(lineSplit[23])
                
                simResults[simId][uid] = [numDescendants, comp, cont]
                
        #print ''
        #print '  Properties of best marker sets:'
        #self.__propertiesOfBestMarkerSets(tree, simResults)
                        
        print(('  Evaluating %d test genomes.' % len(simResults)))
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for testGenomeId in simResults:
            workerQueue.put(testGenomeId)

        for _ in range(numThreads):
            workerQueue.put(None)

        workerProc = [mp.Process(target = self.__workerThread, args = (tree, simResults, metadata, taxonToGenomeIds, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target = self.__writerThread, args = (len(simResults), writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put((None, None, None, None, None, None, None, None, None))
        writeProc.join()
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--threads', help='Threads to use', type=int, default = 40)
    
    args = parser.parse_args()

    markerSetSelection = MarkerSetSelection()
    markerSetSelection.run(args.threads)
