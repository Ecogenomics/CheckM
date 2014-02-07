###############################################################################
#
# treeParser.py - parse genome tree and associated tree metadata 
#
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

import os
import sys
import logging
import json

import dendropy
from  dendropy.dataobject.taxon import Taxon

import prettytable

import defaultValues

from markerSet import MarkerSet

from common import checkDirExists, reassignStdOut, restoreStdOut, getBinIdsFromDir
from seqUtils import readFasta

class TreeParser():
    """Parse genome tree and associated tree metadata."""
    def __init__(self):
        self.logger = logging.getLogger()
        
    def printSummary(self, outputFormat, outDir, resultsParser, bTabTable, outFile):
        if outputFormat == 1:
            self.reportBinTaxonomy(outDir, resultsParser, bTabTable, outFile, bLineageStatistics = False)
        elif outputFormat == 2:
            self.reportBinTaxonomy(outDir, resultsParser, bTabTable, outFile, bLineageStatistics = True)
        elif outputFormat == 3:
            self.reportNewickTree(outDir, outFile, None)
        elif outputFormat == 4:
            self.reportNewickTree(outDir, outFile, 'taxonomy')
        elif outputFormat == 5:
            self.reportFullMSA(outDir, outFile)
        else:
            self.logger.error("Unknown output format: %d", outputFormat)
            
    def reportFullMSA(self, outDir, outFile):
        """Create MSA with all reference and bin alignments."""
    
        # write bin alignments to file
        oldStdOut = reassignStdOut(outFile)
        for line in open(os.path.join(outDir, 'storage', 'tree', defaultValues.PPLACER_CONCAT_SEQ_OUT)):
            print(line.rstrip()) 
        
        # read duplicate seqs
        duplicateNodes = self.__readDuplicateSeqs()
        
        # write reference alignments to file
        seqs = readFasta(os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.concatenated.derep.fasta'))
        for seqId, seq in seqs.iteritems():
            print('>' + seqId)
            print(seq)
                
            if seqId in duplicateNodes:
                for dupSeqId in duplicateNodes[seqId]:
                    print('>' + dupSeqId)
                    print(seq)
                
        restoreStdOut(outFile, oldStdOut)
            
    def readPlacementFile(self, placementFile):
        '''Read pplacer JSON placement file.'''
        jsonData = open(placementFile)
        
        data = json.load(jsonData)
        binIdToPP = {}
        for placementData in data['placements']:
            binId = placementData['nm'][0][0]
            
            topPP = 0
            for pp in placementData['p']:
                if pp[2] > topPP:
                    topPP = pp[2]
                    
            binIdToPP[binId] = topPP
                 
        jsonData.close()
        
        return binIdToPP
            
    def reportNewickTree(self, outDir, outFile, leafLabels=None): 
        # read duplicate nodes
        duplicateSeqs = self.__readDuplicateSeqs()
                         
        # read tree
        treeFile = os.path.join(outDir, 'storage', 'tree', defaultValues.PPLACER_TREE_OUT)  
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # clean up internal node labels
        for node in tree.internal_nodes():
            if node.label:
                labelSplit = node.label.split('|')
                
                label = labelSplit[0]
                if labelSplit[1] != '':
                    label += '|' + labelSplit[1]
                if labelSplit[2] != '':
                    label += '|' + labelSplit[2]
                    
                node.label = label
        
        # insert duplicate nodes into tree
        for leaf in tree.leaf_nodes():
            duplicates = duplicateSeqs.get(leaf.taxon.label, None)
            if duplicates != None:
                newParent = leaf.parent_node.new_child(edge_length = leaf.edge_length)
                curLeaf = leaf.parent_node.remove_child(leaf)
                newParent.new_child(taxon = curLeaf.taxon, edge_length = 0)
                for d in duplicates:
                    newParent.new_child(taxon = Taxon(label = d), edge_length = 0)
                
        # append taxonomy to leaf nodes
        if leafLabels == 'taxonomy':
            # read taxonomy string for each IMG genome
            taxonomy = {}
            for line in open(os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree.taxonomy.tsv')):
                lineSplit = line.split('\t')
                taxonomy[lineSplit[0]] = lineSplit[1].rstrip()
                
            # append taxonomy to leaf labels
            for leaf in tree.leaf_nodes():
                taxaStr = taxonomy.get(leaf.taxon.label, None)
                if taxaStr:
                    leaf.taxon.label += '|' + taxaStr

        # write out tree
        oldStdOut = reassignStdOut(outFile)
        print(tree.as_string(schema='newick', suppress_rooting=True))
        restoreStdOut(outFile, oldStdOut)   
        
    def getBinTaxonomy(self, outDir, binIds): 
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)
               
        # read genome tree
        treeFile = os.path.join(alignOutputDir, defaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # find first parent of each bin with a taxonomic label
        binIdToTaxonomy = {}
        for binId in binIds:
            node = tree.find_node_with_taxon_label(binId)
            if node == None:
                binIdToTaxonomy[binId] = 'NA'
                continue

            # find first node decorated with a taxon string between leaf and root
            taxaStr = None
            parentNode = node.parent_node
            while parentNode != None:     
                if parentNode.label:             
                    tokens = parentNode.label.split('|')
                    
                    if tokens[1] != '':
                        if taxaStr:
                            taxaStr = tokens[1] + ';' + taxaStr
                        else:
                            taxaStr = tokens[1]
                
                parentNode = parentNode.parent_node
            
            if not taxaStr:
                taxaStr = 'k__unclassified'
            binIdToTaxonomy[node.taxon.label] = taxaStr
                
        return binIdToTaxonomy
    
    def getBinMarkerSets(self, outDir, markerFile, 
                                    numGenomesMarkers, numGenomesRefine, 
                                    bootstrap, bNoLineageSpecificRefinement, bRequireTaxonomy):
        """Determine marker set for each bin."""

        self.logger.info('  Determining marker set for each genome bin.')
        
        # get all bin ids
        binIds = getBinIdsFromDir(outDir)
                
        # get statistics for internal nodes
        uniqueIdToLineageStatistics = self.__readNodeMetadata()
                
        # determine marker set for each bin
        treeFile = os.path.join(outDir, 'storage', 'tree', defaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        rootNode = tree.find_node(filter_fn = lambda n: n.parent_node == None)
        
        fout = open(markerFile, 'w')
        fout.write(defaultValues.LINEAGE_MARKER_FILE_HEADER + '\n')
        
        numProcessedBins = 0
        for binId in binIds:
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, len(binIds), float(numProcessedBins)*100/len(binIds))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            # ascend tree to root finding first node meeting all selection criteria
            node = tree.find_node_with_taxon_label(binId)
            if node == None:
                # use universal set
                trustedUniqueId = rootNode.label.split('|')[0]
                trustedStats = uniqueIdToLineageStatistics[trustedUniqueId]
                trustedStats['taxonomy'] = 'root'
            else:
                parentNode = node.parent_node
                while True:
                    if parentNode.label: # nodes inserted by PPLACER will not have a label
                        trustedUniqueId = parentNode.label.split('|')[0]
                        
                        trustedStats = uniqueIdToLineageStatistics[trustedUniqueId]
                        if trustedStats['# genomes'] >= numGenomesMarkers and trustedStats['bootstrap'] >= bootstrap:
                            if not bRequireTaxonomy or trustedStats['taxonomy'] != '':
                                # all criteria meet, so use marker set from this node
                                break
                    
                    if parentNode.parent_node == None:
                        break # reached the root node so terminate
                
                    parentNode = parentNode.parent_node
             
            # get marker set meeting all criteria required for a trusted marker set       
            trustedMarkerSet = eval(trustedStats['marker set'])
            
            # get lineage-specific marker set which will be used to refine the above marker set
            if not bNoLineageSpecificRefinement:
                if node == None:
                    # use universal set
                    uniqueId = rootNode.label.split('|')[0]
                    stats = uniqueIdToLineageStatistics[uniqueId]
                else:
                    parentNode = node.parent_node
                    while True:
                        if parentNode.label: # nodes inserted by PPLACER will not have a label
                            uniqueId = parentNode.label.split('|')[0]
                            stats = uniqueIdToLineageStatistics[uniqueId]
                            
                            if stats['# genomes'] >= numGenomesRefine:
                                break
                        
                        if parentNode.parent_node == None:
                            break # reached the root node so terminate
                    
                        parentNode = parentNode.parent_node
                
                # get lineage-specific marker set       
                lineageMarkerSet = eval(stats['marker set'])
                                        
                # refine marker set by finding the intersection between these two sets,
                # this remove markers that are not single-copy or ubiquitous in the 
                # specific lineage of a bin
                # Note: co-localization information is taken from the trusted set
                
                # get all lineage-specific marker genes
                allLineageSpecificGenes = set()
                for m in lineageMarkerSet:
                    for gene in m:
                        allLineageSpecificGenes.add(gene)
                
                # remove genes not present in the lineage-specific gene set
                finalMarkerSet = []
                for m in trustedMarkerSet:
                    s = set()
                    for gene in m:
                        if gene in allLineageSpecificGenes:
                            s.add(gene)
                                   
                    if s:
                        finalMarkerSet.append(s)
            else:
                # do not refine marker set
                finalMarkerSet = trustedMarkerSet
                
            markerSet = MarkerSet(finalMarkerSet)
            numMarkers, numMarkerSets = markerSet.size()
                
            if bRequireTaxonomy:
                fout.write(binId + '\t' + trustedStats['taxonomy'] + '\t' + str(numMarkers) + '\t' + str(numMarkerSets) + '\t' + str(finalMarkerSet) + '\n')
            else:
                fout.write(binId + '\t' + str(trustedUniqueId) + '\t' + str(numMarkers) + '\t' + str(numMarkerSets) + '\t' + str(finalMarkerSet) + '\n')
                
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
                
        fout.close()
                
    def __readNodeMetadata(self):
        """Read metadata for internal nodes."""
        
        uniqueIdToLineageStatistics = {}
        metadataFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'metadata.tsv')
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
                except:
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
    
    def readLineageMetadata(self, outDir, binIds): 
        """Get metadata for each bin."""
        
        uniqueIdToLineageStatistics = self.__readNodeMetadata()
          
        # read genome tree
        treeFile = os.path.join(outDir, 'storage', 'tree', defaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # find first parent of each bin with a label
        binIdToLineageStatistics = {}
        for binId in binIds:
            node = tree.find_node_with_taxon_label(binId)
            if node == None:
                d = {}
                d['# genomes'] = 'NA'
                d['taxonomy'] = 'NA'          
                d['gc mean'] = 'NA'
                d['gc std'] = 'NA'
                d['genome size mean'] = 'NA'
                d['genome size std'] = 'NA'
                d['gene count mean'] = 'NA'
                d['gene count std'] = 'NA'
                d['marker set'] = 'NA'
                binIdToLineageStatistics[binId] = d
                continue
            
            # find first labeled parent node (nodes inserted by pplacer will be unlabeled)
            parentNode = node.parent_node
            uniqueId = None
            while parentNode != None:
                if parentNode.label:
                    uniqueId = parentNode.label.split('|')[0]
                    break
                    
                parentNode = parentNode.parent_node

            if uniqueId:
                binIdToLineageStatistics[binId] = uniqueIdToLineageStatistics[uniqueId]
            else:
                self.logger.error('Failed to find lineage-specific statistics for inserted bin: ' + node.taxon.label)
                sys.exit(0)
                
        return binIdToLineageStatistics
    
    def reportBinTaxonomy(self, outDir, resultsParser, bTabTable, outFile, bLineageStatistics):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)
        
        # get all bin ids
        binIds = getBinIdsFromDir(outDir)

        # get taxonomy for each bin
        binIdToTaxonomy = self.getBinTaxonomy(outDir, binIds)
        
        # get weighted ML likelihood
        pplacerJsonFile = os.path.join(outDir, 'storage', 'tree', 'concatenated.pplacer.json')
        binIdToWeightedML = self.readPlacementFile(pplacerJsonFile)
        
        # write table
        if not bLineageStatistics:
            self.__printSimpleSummaryTable(binIdToTaxonomy, binIdToWeightedML, resultsParser, bTabTable, outFile)
        else:
            binIdToLineageStatistics = self.readLineageMetadata(outDir, binIds)
            self.__printFullTable(binIdToTaxonomy, binIdToWeightedML, binIdToLineageStatistics, resultsParser, bTabTable, outFile)
        
    def __printSimpleSummaryTable(self, binIdToTaxonomy, binIdToWeightedML, resultsParser, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)

        arbitraryBinId = binIdToTaxonomy.keys()[0]
        header = ['Bin Id', '# marker (of %d)' % len(resultsParser.models[arbitraryBinId]), 'Taxonomy', 'Weighted ML']
        
        if bTabTable: 
            pTable = None
            print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.2'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.align['Taxonomy'] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        for binId in sorted(binIdToTaxonomy.keys()):
            row = [binId, str(len(resultsParser.results[binId].markerHits)), binIdToTaxonomy[binId], binIdToWeightedML.get(binId, 'NA')]
            
            if bTabTable:
                print('\t'.join(map(str, row)))
            else:
                pTable.add_row(row)
                
        if not bTabTable :  
            print(pTable.get_string())
            
        # restore stdout   
        restoreStdOut(outFile, oldStdOut) 
        
    def __printFullTable(self, binIdToTaxonomy, binIdToWeightedML, binIdToLineageStatistics, resultsParser, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)

        arbitraryBinId = binIdToTaxonomy.keys()[0]
        header = ['Bin Id', '# marker (of %d)' % len(resultsParser.models[arbitraryBinId]), 'Taxonomy', 'Weighted ML']
        header += ['# descendant genomes', 'GC mean', 'GC std']
        header += ['Genome size mean (Mbps)', 'Genome size std (Mbps)']
        header += ['Genome count mean', 'Genome count std']
        
        if bTabTable: 
            pTable = None
            print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.2'
            pTable.float_format['GC mean'] = '.1'
            pTable.float_format['GC std'] = '.1'
            pTable.float_format['Genome count mean'] = '.0'
            pTable.float_format['Genome count std'] = '.0'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.align['Taxonomy'] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        for binId in sorted(binIdToTaxonomy.keys()):
            row = [binId, str(len(resultsParser.results[binId].markerHits)), binIdToTaxonomy[binId], binIdToWeightedML.get(binId, 'NA')]
            row += [binIdToLineageStatistics[binId]['# genomes']]
            row += [binIdToLineageStatistics[binId]['gc mean']]
            row += [binIdToLineageStatistics[binId]['gc std']]
            row += [binIdToLineageStatistics[binId]['genome size mean']]
            row += [binIdToLineageStatistics[binId]['genome size std']]
            row += [binIdToLineageStatistics[binId]['gene count mean']]
            row += [binIdToLineageStatistics[binId]['gene count std']]
            
            if bTabTable:
                print('\t'.join(row))
            else:
                pTable.add_row(row)
                
        if not bTabTable :  
            print(pTable.get_string())
            
        # restore stdout   
        restoreStdOut(outFile, oldStdOut)   
        
    def __readDuplicateSeqs(self):
        """Parse file indicating duplicate sequence alignments."""
        duplicateSeqs = {}
        for line in open(os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree.derep.txt')):
            lineSplit = line.rstrip().split()
            if len(lineSplit) > 1:
                duplicateSeqs[lineSplit[0]] = lineSplit[1:]
                
        return duplicateSeqs
    