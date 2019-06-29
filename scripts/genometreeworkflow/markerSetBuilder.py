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
Methods for identifying and investigating marker sets.
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
import random
from collections import defaultdict
from checkm.markerSets import BinMarkerSets, MarkerSet

from checkm.util.img import IMG

from numpy import array
from numpy.random import choice

class MarkerSetBuilder(object):
    def __init__(self):
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        self.colocatedFile = './data/colocated.tsv'
        self.duplicateSeqs = self.readDuplicateSeqs()
        self.uniqueIdToLineageStatistics = self.__readNodeMetadata()
        
        self.cachedGeneCountTable = None
        
    def precomputeGenomeSeqLens(self, genomeIds):
        """Cache the length of contigs/scaffolds for all genomes."""
        
        # This function is intended to speed up functions, such as img.geneDistTable(),
        # that are called multiple times (typically during simulations)
        self.img.precomputeGenomeSeqLens(genomeIds)
            
    def precomputeGenomeFamilyPositions(self, genomeIds, spacingBetweenContigs):
        """Cache position of PFAM and TIGRFAM genes in genomes."""
        
        # This function is intended to speed up functions, such as img.geneDistTable(),
        # that are called multiple times (typically during simulations)
        self.img.precomputeGenomeFamilyPositions(genomeIds, spacingBetweenContigs)
        
    def precomputeGenomeFamilyScaffolds(self, genomeIds):
        """Cache scaffolds of PFAM and TIGRFAM genes in genomes."""
        
        # This function is intended to speed up functions, such as img.geneDistTable(),
        # that are called multiple times (typically during simulations)
        self.img.precomputeGenomeFamilyScaffolds(genomeIds)
        
    def getLineageMarkerGenes(self, lineage, minGenomes = 20, minMarkerSets = 20):
        pfamIds = set()
        tigrIds = set()

        bHeader = True
        for line in open(self.colocatedFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            curLineage = lineSplit[0]
            numGenomes = int(lineSplit[1])
            numMarkerSets = int(lineSplit[3])
            markerSets = lineSplit[4:]

            if curLineage != lineage or numGenomes < minGenomes or numMarkerSets < minMarkerSets:
                continue

            for ms in markerSets:
                markers = ms.split(',')
                for m in markers:
                    if 'pfam' in m:
                        pfamIds.add(m.strip())
                    elif 'TIGR' in m:
                        tigrIds.add(m.strip())

        return pfamIds, tigrIds

    def getCalculatedMarkerGenes(self, minGenomes = 20, minMarkerSets = 20):
        pfamIds = set()
        tigrIds = set()

        bHeader = True
        for line in open(self.colocatedFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            numGenomes = int(lineSplit[1])
            numMarkerSets = int(lineSplit[3])
            markerSets = lineSplit[4:]

            if numGenomes < minGenomes or numMarkerSets < minMarkerSets:
                continue

            for ms in markerSets:
                markers = ms.split(',')
                for m in markers:
                    if 'pfam' in m:
                        pfamIds.add(m.strip())
                    elif 'TIGR' in m:
                        tigrIds.add(m.strip())

        return pfamIds, tigrIds

    def markerGenes(self, genomeIds, countTable, ubiquityThreshold, singleCopyThreshold):
        if ubiquityThreshold < 1 or singleCopyThreshold < 1:
            print '[Warning] Looks like degenerate threshold.'

        # find genes meeting ubiquity and single-copy thresholds
        markers = set()
        for clusterId, genomeCounts in countTable.iteritems():
            ubiquity = 0
            singleCopy = 0
            
            if len(genomeCounts) < ubiquityThreshold:
                # gene is clearly not ubiquitous
                continue
            
            for genomeId in genomeIds:
                count = genomeCounts.get(genomeId, 0)

                if count > 0:
                    ubiquity += 1

                if count == 1:
                    singleCopy += 1

            if ubiquity >= ubiquityThreshold and singleCopy >= singleCopyThreshold:
                markers.add(clusterId)

        return markers

    def colocatedGenes(self, geneDistTable, distThreshold = 5000, genomeThreshold = 0.95):
        """Identify co-located gene pairs."""
                
        colocatedGenes = defaultdict(int)
        for _, clusterIdToGeneLocs in geneDistTable.iteritems():
            clusterIds = clusterIdToGeneLocs.keys()
            for i, clusterId1 in enumerate(clusterIds):
                geneLocations1 = clusterIdToGeneLocs[clusterId1]
                
                for clusterId2 in clusterIds[i+1:]:
                    geneLocations2 = clusterIdToGeneLocs[clusterId2]
                    bColocated = False
                    for p1 in geneLocations1:
                        for p2 in geneLocations2:
                            if abs(p1[0] - p2[0]) < distThreshold:
                                bColocated = True
                                break

                        if bColocated:
                            break

                    if bColocated:
                        if clusterId1 <= clusterId2:
                            colocatedStr = clusterId1 + '-' + clusterId2
                        else:
                            colocatedStr = clusterId2 + '-' + clusterId1
                        colocatedGenes[colocatedStr] += 1

        colocated = []
        for colocatedStr, count in colocatedGenes.iteritems():
            if float(count)/len(geneDistTable) > genomeThreshold:
                colocated.append(colocatedStr)

        return colocated

    def colocatedSets(self, colocatedGenes, markerGenes):
        # run through co-located genes once creating initial sets
        sets = []
        for cg in colocatedGenes:
            geneA, geneB = cg.split('-')
            sets.append(set([geneA, geneB]))

        # combine any sets with overlapping genes
        bProcessed = [False]*len(sets)
        finalSets = []
        for i in xrange(0, len(sets)):
            if bProcessed[i]:
                continue

            curSet = sets[i]
            bProcessed[i] = True

            bUpdated = True
            while bUpdated:
                bUpdated = False
                for j in xrange(i+1, len(sets)):
                    if bProcessed[j]:
                        continue

                    if len(curSet.intersection(sets[j])) > 0:
                        curSet.update(sets[j])
                        bProcessed[j] = True
                        bUpdated = True

            finalSets.append(curSet)

        # add all singletons into colocated sets
        for clusterId in markerGenes:
            bFound = False
            for cs in finalSets:
                if clusterId in cs:
                    bFound = True

            if not bFound:
                finalSets.append(set([clusterId]))

        return finalSets

    def genomeCheck(self, colocatedSet, genomeId, countTable):
        comp = 0.0
        cont = 0.0
        missingMarkers = set()
        duplicateMarkers = set()

        if len(colocatedSet) == 0:
            return comp, cont, missingMarkers, duplicateMarkers

        for cs in colocatedSet:
            present = 0
            multiCopy = 0
            for contigId in cs:
                count = countTable[contigId].get(genomeId, 0)
                if count == 1:
                    present += 1
                elif count > 1:
                    present += 1
                    multiCopy += (count-1)
                    duplicateMarkers.add(contigId)
                elif count == 0:
                    missingMarkers.add(contigId)

            comp += float(present) / len(cs)
            cont += float(multiCopy) / len(cs)

        return comp / len(colocatedSet), cont / len(colocatedSet), missingMarkers, duplicateMarkers

    def uniformity(self, genomeSize, pts):
        U = float(genomeSize) / (len(pts)+1)  # distance between perfectly evenly spaced points

        # calculate distance between adjacent points
        dists = []
        pts = sorted(pts)
        for i in xrange(0, len(pts)-1):
            dists.append(pts[i+1] - pts[i])

        # calculate uniformity index
        num = 0
        den = 0
        for d in dists:
            num += abs(d - U)
            den += max(d, U)

        return 1.0 - num/den
    
    def sampleGenome(self, genomeLen, percentComp, percentCont, contigLen):
        """Sample a genome to simulate a given percent completion and contamination."""
        
        contigsInGenome = genomeLen / contigLen

        # determine number of contigs to achieve desired completeness and contamination
        contigsToSampleComp = int(contigsInGenome*percentComp + 0.5)
        contigsToSampleCont = int(contigsInGenome*percentCont + 0.5)

        # randomly sample contigs with contamination done via sampling with replacement
        compContigs = random.sample(xrange(contigsInGenome), contigsToSampleComp)  
        contContigs = choice(xrange(contigsInGenome), contigsToSampleCont, replace=True)
    
        # determine start of each contig
        contigStarts = [c*contigLen for c in compContigs]
        contigStarts += [c*contigLen for c in contContigs]
            
        contigStarts.sort()
        
        trueComp = float(contigsToSampleComp)*contigLen*100 / genomeLen
        trueCont = float(contigsToSampleCont)*contigLen*100 / genomeLen

        return trueComp, trueCont, contigStarts
    
    def sampleGenomeScaffoldsInvLength(self, targetPer, seqLens, genomeSize):
        """Sample genome comprised of several sequences with probability inversely proportional to length."""
        
        # calculate probability of sampling a sequences
        seqProb = []
        for _, seqLen in seqLens.iteritems():
            prob = 1.0 / (float(seqLen) / genomeSize)
            seqProb.append(prob)
            
        seqProb = array(seqProb)
        seqProb /= sum(seqProb)
            
        # select sequence with probability proportional to length
        selectedSeqsIds = choice(seqLens.keys(), size = len(seqLens), replace=False, p = seqProb)
        
        sampledSeqIds = []
        truePer = 0.0
        for seqId in selectedSeqsIds:
            sampledSeqIds.append(seqId)
            truePer += float(seqLens[seqId]) / genomeSize
            
            if truePer >= targetPer:
                break
        
        return sampledSeqIds, truePer*100
    
    def sampleGenomeScaffoldsWithoutReplacement(self, targetPer, seqLens, genomeSize):
        """Sample genome comprised of several sequences without replacement.
        
          Sampling is conducted randomly until the selected sequences comprise
          greater than or equal to the desired target percentage.
        """
 
        selectedSeqsIds = choice(seqLens.keys(), size = len(seqLens), replace=False)
        
        sampledSeqIds = []
        truePer = 0.0
        for seqId in selectedSeqsIds:
            sampledSeqIds.append(seqId)
            truePer += float(seqLens[seqId]) / genomeSize
            
            if truePer >= targetPer:
                break
        
        return sampledSeqIds, truePer*100
    
    def containedMarkerGenes(self, markerGenes, clusterIdToGenomePositions, startPartialGenomeContigs, contigLen):
        """Determine markers contained in a set of contigs."""
        
        contained = {}
        for markerGene in markerGenes:
            positions = clusterIdToGenomePositions.get(markerGene, [])

            containedPos = []
            for p in positions:
                for s in startPartialGenomeContigs:
                    if (p[0] - s) >= 0 and (p[0] - s) < contigLen:
                        containedPos.append(s)

            if len(containedPos) > 0:
                contained[markerGene] = containedPos

        return contained
    
    def markerGenesOnScaffolds(self, markerGenes, genomeId, scaffoldIds, containedMarkerGenes):
        """Determine if marker genes are found on the scaffolds of a given genome."""
        for markerGeneId in markerGenes:
            scaffoldIdsWithMarker = self.img.cachedGenomeFamilyScaffolds[genomeId].get(markerGeneId, [])

            for scaffoldId in scaffoldIdsWithMarker:
                if scaffoldId in scaffoldIds:
                    containedMarkerGenes[markerGeneId] += [scaffoldId]
        
    def readDuplicateSeqs(self):
        """Parse file indicating duplicate sequence alignments."""
        duplicateSeqs = {}
        for line in open(os.path.join('/srv/whitlam/bio/db/checkm/genome_tree', 'genome_tree.derep.txt')):
            lineSplit = line.rstrip().split()
            if len(lineSplit) > 1:
                duplicateSeqs[lineSplit[0]] = lineSplit[1:]
                
        return duplicateSeqs
        
    def __readNodeMetadata(self):
        """Read metadata for internal nodes."""
        
        uniqueIdToLineageStatistics = {}
        metadataFile = os.path.join('/srv/whitlam/bio/db/checkm/genome_tree', 'genome_tree.metadata.tsv')
        with open(metadataFile) as f:
            f.readline()
            for line in f:
                lineSplit = line.rstrip().split('\t')
                
                uniqueId = lineSplit[0]
                
                d = {}
                d['# genomes'] = int(lineSplit[1])
                d['taxonomy'] = lineSplit[2]
                try:
                    d['bootstrap'] = float(lineSplit[3])
                except ValueError:
                    d['bootstrap'] = 'NA'                 
                d['gc mean'] = float(lineSplit[4])
                d['gc std'] = float(lineSplit[5])
                d['genome size mean'] = float(lineSplit[6])/1e6
                d['genome size std'] = float(lineSplit[7])/1e6
                d['gene count mean'] = float(lineSplit[8])
                d['gene count std'] = float(lineSplit[9])
                d['marker set'] = lineSplit[10].rstrip()
                
                uniqueIdToLineageStatistics[uniqueId] = d
                
        return uniqueIdToLineageStatistics
    
    def __getNextNamedNode(self, node, uniqueIdToLineageStatistics):
        """Get first parent node with taxonomy information."""
        parentNode = node.parent_node
        while True:
            if parentNode == None:
                break # reached the root node so terminate
            
            if parentNode.label:
                trustedUniqueId = parentNode.label.split('|')[0]
                trustedStats = uniqueIdToLineageStatistics[trustedUniqueId]
                if trustedStats['taxonomy'] != '':
                    return trustedStats['taxonomy']
                
            parentNode = parentNode.parent_node            
                    
        return 'root'
    
    def __refineMarkerSet(self, markerSet, lineageSpecificMarkerSet):
        """Refine marker set to account for lineage-specific gene loss and duplication."""
                                        
        # refine marker set by finding the intersection between these two sets,
        # this removes markers that are not single-copy or ubiquitous in the 
        # specific lineage of a bin
        # Note: co-localization information is taken from the trusted set
                
        # remove genes not present in the lineage-specific gene set
        finalMarkerSet = []
        for ms in markerSet.markerSet:
            s = set()
            for gene in ms:
                if gene in lineageSpecificMarkerSet.getMarkerGenes():
                    s.add(gene)
                           
            if s:
                finalMarkerSet.append(s)

        refinedMarkerSet = MarkerSet(markerSet.UID, markerSet.lineageStr, markerSet.numGenomes, finalMarkerSet)
    
        return refinedMarkerSet
    
    def ____removeInvalidLineageMarkerGenes(self, markerSet, lineageSpecificMarkersToRemove):
        """Refine marker set to account for lineage-specific gene loss and duplication."""
                                        
        # refine marker set by removing marker genes subject to lineage-specific
        # gene loss and duplication 
        #
        # Note: co-localization information is taken from the trusted set
                
        finalMarkerSet = []
        for ms in markerSet.markerSet:
            s = set()
            for gene in ms:
                if gene.startswith('PF'):
                    print 'ERROR! Expected genes to start with pfam, not PF.'
                    
                if gene not in lineageSpecificMarkersToRemove:
                    s.add(gene)
                           
            if s:
                finalMarkerSet.append(s)

        refinedMarkerSet = MarkerSet(markerSet.UID, markerSet.lineageStr, markerSet.numGenomes, finalMarkerSet)
    
        return refinedMarkerSet
    
    def missingGenes(self, genomeIds, markerGenes, ubiquityThreshold):
        """Inferring consistently missing marker genes within a set of genomes."""
        
        if self.cachedGeneCountTable != None:
            geneCountTable = self.cachedGeneCountTable
        else:
            geneCountTable = self.img.geneCountTable(genomeIds)
        
        # find genes meeting ubiquity and single-copy thresholds
        missing = set()
        for clusterId, genomeCounts in geneCountTable.iteritems():
            if clusterId not in markerGenes:
                continue
                 
            absence = 0 
            for genomeId in genomeIds:
                count = genomeCounts.get(genomeId, 0)

                if count == 0:
                    absence += 1

            if absence >= ubiquityThreshold*len(genomeIds):
                missing.add(clusterId)

        return missing
    
    def duplicateGenes(self, genomeIds, markerGenes, ubiquityThreshold):
        """Inferring consistently duplicated marker genes within a set of genomes."""
        
        if self.cachedGeneCountTable != None:
            geneCountTable = self.cachedGeneCountTable
        else:
            geneCountTable = self.img.geneCountTable(genomeIds)
        
        # find genes meeting ubiquity and single-copy thresholds
        duplicate = set()
        for clusterId, genomeCounts in geneCountTable.iteritems():
            if clusterId not in markerGenes:
                continue
                 
            duplicateCount = 0 
            for genomeId in genomeIds:
                count = genomeCounts.get(genomeId, 0)

                if count > 1:
                    duplicateCount += 1

            if duplicateCount >= ubiquityThreshold*len(genomeIds):
                duplicate.add(clusterId)

        return duplicate
    
    def buildMarkerGenes(self, genomeIds, ubiquityThreshold, singleCopyThreshold):
        """Infer marker genes from specified genomes."""
        
        if self.cachedGeneCountTable != None:
            geneCountTable = self.cachedGeneCountTable
        else:
            geneCountTable = self.img.geneCountTable(genomeIds)
        
        #counts = []
        #singleCopy = 0
        #for genomeId, count in geneCountTable['pfam01351'].iteritems():
        #    print genomeId, count
        #    counts.append(count)
        #    if count == 1:
        #        singleCopy += 1
            
        #print 'Ubiquity: %d of %d' % (len(counts), len(genomeIds))
        #print 'Single-copy: %d of %d' % (singleCopy, len(genomeIds))
        #print 'Mean: %.2f' % mean(counts)

        markerGenes = self.markerGenes(genomeIds, geneCountTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
        tigrToRemove = self.img.identifyRedundantTIGRFAMs(markerGenes)
        markerGenes = markerGenes - tigrToRemove

        return markerGenes
    
    def buildMarkerSet(self, genomeIds, ubiquityThreshold, singleCopyThreshold, spacingBetweenContigs = 5000):  
        """Infer marker set from specified genomes."""      

        markerGenes = self.buildMarkerGenes(genomeIds, ubiquityThreshold, singleCopyThreshold)

        geneDistTable = self.img.geneDistTable(genomeIds, markerGenes, spacingBetweenContigs)
        colocatedGenes = self.colocatedGenes(geneDistTable)
        colocatedSets = self.colocatedSets(colocatedGenes, markerGenes)
        markerSet = MarkerSet(0, 'NA', len(genomeIds), colocatedSets)

        return markerSet
    
    def readLineageSpecificGenesToRemove(self):
        """Get set of genes subject to lineage-specific gene loss and duplication."""       
    
        self.lineageSpecificGenesToRemove = {}
        for line in open('/srv/whitlam/bio/db/checkm/genome_tree/missing_duplicate_genes_50.tsv'):
            lineSplit = line.split('\t')
            uid = lineSplit[0]
            missingGenes = eval(lineSplit[1])
            duplicateGenes = eval(lineSplit[2])
            self.lineageSpecificGenesToRemove[uid] = missingGenes.union(duplicateGenes)
            
    def buildBinMarkerSet(self, tree, curNode, ubiquityThreshold, singleCopyThreshold, bMarkerSet = True, genomeIdsToRemove = None):   
        """Build lineage-specific marker sets for a genome in a LOO-fashion."""
                               
        # determine marker sets for bin      
        binMarkerSets = BinMarkerSets(curNode.label, BinMarkerSets.TREE_MARKER_SET)
        refinedBinMarkerSet = BinMarkerSets(curNode.label, BinMarkerSets.TREE_MARKER_SET)         

        # ascend tree to root, recording all marker sets 
        uniqueId = curNode.label.split('|')[0] 
        lineageSpecificRefinement = self.lineageSpecificGenesToRemove[uniqueId]
        
        while curNode != None:
            uniqueId = curNode.label.split('|')[0] 
            stats = self.uniqueIdToLineageStatistics[uniqueId]
            taxonomyStr = stats['taxonomy']
            if taxonomyStr == '':
                taxonomyStr = self.__getNextNamedNode(curNode, self.uniqueIdToLineageStatistics)

            leafNodes = curNode.leaf_nodes()
            genomeIds = set()
            for leaf in leafNodes:
                genomeIds.add(leaf.taxon.label.replace('IMG_', ''))
                
                duplicateGenomes = self.duplicateSeqs.get(leaf.taxon.label, [])
                for dup in duplicateGenomes:
                    genomeIds.add(dup.replace('IMG_', ''))

            # remove all genomes from the same taxonomic group as the genome of interest
            if genomeIdsToRemove != None:
                genomeIds.difference_update(genomeIdsToRemove) 

            if len(genomeIds) >= 2:
                if bMarkerSet:
                    markerSet = self.buildMarkerSet(genomeIds, ubiquityThreshold, singleCopyThreshold)
                else:
                    markerSet = MarkerSet(0, 'NA', len(genomeIds), [self.buildMarkerGenes(genomeIds, ubiquityThreshold, singleCopyThreshold)])
                
                markerSet.lineageStr = uniqueId + ' | ' + taxonomyStr.split(';')[-1]
                binMarkerSets.addMarkerSet(markerSet)
        
                #refinedMarkerSet = self.__refineMarkerSet(markerSet, lineageSpecificMarkerSet)
                refinedMarkerSet = self.____removeInvalidLineageMarkerGenes(markerSet, lineageSpecificRefinement)
                #print 'Refinement: %d of %d' % (len(refinedMarkerSet.getMarkerGenes()), len(markerSet.getMarkerGenes()))
                refinedBinMarkerSet.addMarkerSet(refinedMarkerSet)
            
            curNode = curNode.parent_node
                
        return binMarkerSets, refinedBinMarkerSet
    
    def buildDomainMarkerSet(self, tree, curNode, ubiquityThreshold, singleCopyThreshold, bMarkerSet = True, genomeIdsToRemove = None):   
        """Build domain-specific marker sets for a genome in a LOO-fashion."""
                               
        # determine marker sets for bin      
        binMarkerSets = BinMarkerSets(curNode.label, BinMarkerSets.TREE_MARKER_SET)
        refinedBinMarkerSet = BinMarkerSets(curNode.label, BinMarkerSets.TREE_MARKER_SET)         

        # calculate marker set for bacterial or archaeal node
        uniqueId = curNode.label.split('|')[0] 
        lineageSpecificRefinement = self.lineageSpecificGenesToRemove[uniqueId]
        
        while curNode != None:
            uniqueId = curNode.label.split('|')[0] 
            if uniqueId != 'UID2' and uniqueId != 'UID203':
                curNode = curNode.parent_node
                continue

            stats = self.uniqueIdToLineageStatistics[uniqueId]
            taxonomyStr = stats['taxonomy']
            if taxonomyStr == '':
                taxonomyStr = self.__getNextNamedNode(curNode, self.uniqueIdToLineageStatistics)

            leafNodes = curNode.leaf_nodes()
            genomeIds = set()
            for leaf in leafNodes:
                genomeIds.add(leaf.taxon.label.replace('IMG_', ''))
                
                duplicateGenomes = self.duplicateSeqs.get(leaf.taxon.label, [])
                for dup in duplicateGenomes:
                    genomeIds.add(dup.replace('IMG_', ''))

            # remove all genomes from the same taxonomic group as the genome of interest
            if genomeIdsToRemove != None:
                genomeIds.difference_update(genomeIdsToRemove) 

            if len(genomeIds) >= 2:
                if bMarkerSet:
                    markerSet = self.buildMarkerSet(genomeIds, ubiquityThreshold, singleCopyThreshold)
                else:
                    markerSet = MarkerSet(0, 'NA', len(genomeIds), [self.buildMarkerGenes(genomeIds, ubiquityThreshold, singleCopyThreshold)])
                
                markerSet.lineageStr = uniqueId + ' | ' + taxonomyStr.split(';')[-1]
                binMarkerSets.addMarkerSet(markerSet)
        
                #refinedMarkerSet = self.__refineMarkerSet(markerSet, lineageSpecificMarkerSet)
                refinedMarkerSet = self.____removeInvalidLineageMarkerGenes(markerSet, lineageSpecificRefinement)
                #print 'Refinement: %d of %d' % (len(refinedMarkerSet.getMarkerGenes()), len(markerSet.getMarkerGenes()))
                refinedBinMarkerSet.addMarkerSet(refinedMarkerSet)
            
            curNode = curNode.parent_node
                
        return binMarkerSets, refinedBinMarkerSet