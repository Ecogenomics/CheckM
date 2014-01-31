###############################################################################
#
# pplacer.py - runs pplacer and provides functions for parsing output  
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
import subprocess
import logging
import json
from collections import defaultdict

import dendropy
from  dendropy.dataobject.taxon import Taxon

import prettytable

import defaultValues

from common import checkDirExists, reassignStdOut, restoreStdOut
from seqUtils import readFasta, writeFasta

class PplacerRunner():
    """Wrapper for running pplacer."""
    def __init__(self, threads):
        self.logger = logging.getLogger()
        self.numThreads = threads
        
        # make sure pplace and guppy are on the system path
        self.__checkForPplacer()
        self.__checkForGuppy()
        
        self.CONCAT_SEQ_OUT = 'concatenated.fasta'
        self.PPLACER_JSON_OUT = 'concatenated.pplacer.json'
        self.PPLACER_OUT = 'pplacer.out'
        self.TREE_OUT = 'concatenated.tre'
          
    def run(self, binFiles, outDir):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        # create concatenated alignment file for each bin
        concatenatedAlignFile = self.__createConcatenatedAlignment(binFiles, alignOutputDir)
        
        # run pplacer to place bins in reference genome tree
        self.logger.info('  Placing %d bins in the genome tree with pplacer.' % len(binFiles))
        refpkg = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree_prok.refpkg')
        pplacerJsonOut = os.path.join(alignOutputDir, self.PPLACER_JSON_OUT)
        pplacerOut = os.path.join(alignOutputDir, self.PPLACER_OUT)
        cmd = 'pplacer -j %d -c %s -o %s %s > %s' % (self.numThreads, refpkg, pplacerJsonOut, concatenatedAlignFile, pplacerOut)
        os.system(cmd)
        
        # extract tree
        treeFile = os.path.join(alignOutputDir, self.TREE_OUT)
        cmd = 'guppy tog -o %s %s' % (treeFile, pplacerJsonOut)
        os.system(cmd)
        
    def printSummary(self, outputFormat, outDir, resultsParser, bTabTable, outFile):
        if outputFormat == 1:
            self.reportBinTaxonomy(outDir, resultsParser, bTabTable, outFile)
        elif outputFormat == 2:
            self.reportNewickTree(outDir, outFile, None)
        elif outputFormat == 3:
            self.reportNewickTree(outDir, outFile, 'taxonomy')
        elif outputFormat == 4:
            self.reportFullMSA(outDir, outFile)
        else:
            self.logger.error("Unknown output format: %d", outputFormat)
            
    def reportFullMSA(self, outDir, outFile):
        """Create MSA with all reference and bin alignments."""
    
        # write bin alignments to file
        oldStdOut = reassignStdOut(outFile)
        for line in open(os.path.join(outDir, 'storage', 'tree', self.CONCAT_SEQ_OUT)):
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
        treeFile = os.path.join(outDir, 'storage', 'tree', self.TREE_OUT)  
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
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
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        treeFile = os.path.join(alignOutputDir, self.TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # find first parent of each bin with a taxonomic label
        binIdToTaxonomy = {}
        for node in tree.leaf_nodes():
            if node.taxon.label in binIds:
                # find first node decorated with a taxon string between leaf and root
                taxaStr = None
                parentNode = node.parent_node
                while parentNode != None:                  
                    if parentNode.label and '|' in parentNode.label:
                        tokens = parentNode.label.split('|')
                        
                        if taxaStr:
                            taxaStr = tokens[0] + ';' + taxaStr
                        else:
                            taxaStr = tokens[0]
                    
                    parentNode = parentNode.parent_node
                    if parentNode == None:
                        break
                
                if not taxaStr:
                    taxaStr = 'k__unclassified'
                binIdToTaxonomy[node.taxon.label] = taxaStr
                
        return binIdToTaxonomy
    
    def reportBinTaxonomy(self, outDir, resultsParser, bTabTable, outFile):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)
        
        # get all bin ids
        files = os.listdir(outDir)
        
        binIds = set()
        for f in files:
            if os.path.isdir(os.path.join(outDir, f)) and f != 'storage':
                binIds.add(f)
        
        # get taxonomy for each bin
        binIdToTaxonomy = self.getBinTaxonomy(outDir, binIds)
        
        # get weighted ML likelihood
        pplacerJsonFile = os.path.join(outDir, 'storage', 'tree', 'concatenated.pplacer.json')
        binIdToWeightedML = self.readPlacementFile(pplacerJsonFile)
        
        # write table
        self.__printTable(binIdToTaxonomy, binIdToWeightedML, resultsParser, bTabTable, outFile)
        
    def __printTable(self, binIdToTaxonomy, binIdToWeightedML, resultsParser, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)

        header = ['Bin Id', '# marker (of %d)' % len(resultsParser.models), 'Taxonomy', 'Weighted ML']
        
        if bTabTable: 
            pTable = None
            print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.3'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.align['Taxonomy'] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        for binId in sorted(binIdToTaxonomy.keys()):
            row = [binId, str(len(resultsParser.results[binId].markerHits)), binIdToTaxonomy[binId], binIdToWeightedML[binId]]
            
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

    def __createConcatenatedAlignment(self, binFiles, alignOutputDir):
        """Create a concatenated alignment of marker genes for each bin."""
                
        # read alignment files
        self.logger.info('  Reading marker alignment files.')
        alignments = defaultdict(dict)
        files = os.listdir(alignOutputDir)
        binIds = set()
        for f in files:
            if f.endswith('.masked.faa'):
                markerId = f[0:f.find('.')]
                seqs = readFasta(os.path.join(alignOutputDir, f))
                
                for seqId, seq in seqs.iteritems():
                    binId = seqId[0:seqId.find(defaultValues.SEQ_CONCAT_CHAR)]

                    alignments[markerId][binId] = seq
                    binIds.add(binId)

        # create concatenated alignment
        self.logger.info('  Concatenating alignments.')
        concatenatedSeqs = {}
        t = 0
        for markerId in sorted(alignments.keys()):
            seqs = alignments[markerId]
            alignLen = len(seqs[seqs.keys()[0]])
            
            t += alignLen

            for binId in binIds:
                if binId in seqs:
                    # append alignment
                    concatenatedSeqs[binId] = concatenatedSeqs.get(binId, '') + seqs[binId]
                else:
                    # missing gene
                    concatenatedSeqs[binId] = concatenatedSeqs.get(binId, '') + '-'*alignLen

        # save concatenated alignment
        concatenatedAlignFile = os.path.join(alignOutputDir, self.CONCAT_SEQ_OUT)
        writeFasta(concatenatedSeqs, concatenatedAlignFile)
        
        return concatenatedAlignFile

    def __checkForPplacer(self):
        """Check to see if pplacer is on the system before we try to run it."""
        
        # Assume that a successful pplacer -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['pplacer', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure pplacer is on your system path.")
            sys.exit()
            
    def __checkForGuppy(self):
        """Check to see if guppy is on the system before we try to run it."""
        
        # Assume that a successful pplacer -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['guppy', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure guppy is on your system path.")
            sys.exit()

class PplacerParser():
    """Parses pplacer output."""
    def __init__(self):
        pass

    
