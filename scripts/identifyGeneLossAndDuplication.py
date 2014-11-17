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
Identify marker genes that are lost or duplicated within a lineage.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import time
import argparse

from checkm.util.img import IMG
from lib.markerSetBuilder import MarkerSetBuilder

from checkm.treeParser import TreeParser
from checkm.markerSets import MarkerSet

import dendropy

class IdentifyGeneLossAndDuplication(object):
    def __init__(self):
        self.markerSetBuilder = MarkerSetBuilder()
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')

    def run(self, ubiquityThreshold, minGenomes):
        # Pre-compute gene count table
        print 'Computing gene count table.'
        start = time.time()
        metadata = self.img.genomeMetadata()
        self.markerSetBuilder.cachedGeneCountTable = self.img.geneCountTable(metadata.keys())
        end = time.time()
        print '    globalGeneCountTable: %.2f' % (end - start)

        # read selected node for defining marker set
        print 'Reading node defining marker set for each internal node.'
        selectedMarkerNode = {}
        for line in open('/srv/whitlam/bio/db/checkm/selected_marker_sets.tsv'):
            lineSplit = line.split('\t')
            selectedMarkerNode[lineSplit[0].strip()] = lineSplit[1].strip()
            
        # read duplicate taxa
        print 'Reading list of identical taxa in genome tree.'
        duplicateTaxa = {}
        for line in open('/srv/whitlam/bio/db/checkm/genome_tree/genome_tree.derep.txt'):
            lineSplit = line.rstrip().split()
            if len(lineSplit) > 1:
                duplicateTaxa[lineSplit[0]] = lineSplit[1:]
        
        # read in node metadata
        print 'Reading node metadata.'
        treeParser = TreeParser()
        uniqueIdToLineageStatistics = treeParser.readNodeMetadata()
        
        # read genome tree
        print 'Reading in genome tree.'
                
        treeFile = '/srv/whitlam/bio/db/checkm/genome_tree/genome_tree_prok.refpkg/genome_tree.final.tre'
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # determine lineage-specific gene loss and duplication (relative to potential marker genes used by a node)
        print 'Determining lineage-specific gene loss and duplication'
        
        fout = open('/srv/whitlam/bio/db/checkm/genome_tree/missing_duplicate_genes_50.tsv', 'w')
        
        processed = 0
        numInternalNodes = len(tree.internal_nodes())
        for node in tree.internal_nodes():
            processed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) internal nodes.' % (processed, numInternalNodes, float(processed)*100/numInternalNodes)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            nodeId = node.label.split('|')[0]
            
            missingGenes = []
            duplicateGenes = []
            
            nodeStats = uniqueIdToLineageStatistics[nodeId]
            if nodeStats['# genomes'] >= minGenomes:               
                # get marker genes defined for current node along with all parental nodes    
                markerGenes = set() 
                parentNode = node
                while parentNode != None:                     
                    parentNodeId = parentNode.label.split('|')[0]
                
                    stats = uniqueIdToLineageStatistics[parentNodeId]
                    markerSet = MarkerSet(parentNodeId, stats['taxonomy'], stats['# genomes'], eval(stats['marker set']))
                    markerGenes = markerGenes.union(markerSet.getMarkerGenes())
                
                    parentNode = parentNode.parent_node
                
                # silly hack since PFAM ids are inconsistent between the PFAM data and IMG data
                revisedMarkerGeneIds = set()
                for mg in markerGenes:
                    if mg.startswith('PF'):
                        revisedMarkerGeneIds.add(mg[0:mg.rfind('.')].replace('PF', 'pfam'))
                    else:
                        revisedMarkerGeneIds.add(mg)
                
                # get all genomes below the internal node (including genomes removed as duplicates)
                genomeIds = []
                for leaf in node.leaf_nodes():
                    genomeIds.append(leaf.taxon.label.replace('IMG_', ''))
                    if leaf.taxon.label in duplicateTaxa:
                        for genomeId in duplicateTaxa[leaf.taxon.label]:
                            genomeIds.append(genomeId.replace('IMG_', ''))
                            
                    genomeIds.append(leaf.taxon.label.replace('IMG_', ''))
                
                missingGenes = self.markerSetBuilder.missingGenes(genomeIds, revisedMarkerGeneIds, ubiquityThreshold)
                duplicateGenes = self.markerSetBuilder.duplicateGenes(genomeIds, revisedMarkerGeneIds, ubiquityThreshold)
                
            fout.write('%s\t%s\t%s\n' % (nodeId, str(missingGenes), str(duplicateGenes)))
            
        sys.stdout.write('\n')
            
        fout.close()
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify co-located genes within genomes for a specific lineage.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='ubiquity threshold for defining gene loss or duplication', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='minimum genomes required to include in analysis', type=int, default = 2)

    args = parser.parse_args()

    identifyGeneLossAndDuplication = IdentifyGeneLossAndDuplication()
    identifyGeneLossAndDuplication.run(args.ubiquity, args.min_genomes)
