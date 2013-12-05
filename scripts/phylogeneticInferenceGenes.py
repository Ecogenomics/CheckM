#!/usr/bin/env python

###############################################################################
#
# phylogeneticInferenceGenes - identify genes suitable for phylogenetic inference
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

from lib.img import IMG
from lib.markerSet import MarkerSet

from checkm.seqUtils import readFasta
from checkm.hmmer import HMMERRunner

class PhylogeneticInferenceGenes(object):
    def __init__(self):
        pass
    
    def run(self, phyloUbiquityThreshold, phyloSingleCopyThreshold, numThreads):   
        img = IMG()
        markerset = MarkerSet()
        
        metadata = img.genomeMetadata()
          
        allTrustedGenomeIds = set()
        phyloMarkerGenes = {}
        for lineage in ['Archaea', 'Bacteria']:  
            # get all genomes in lineage
            print '\nIdentifying all ' + lineage + ' genomes.'
            trustedGenomeIds = img.genomeIdsByTaxonomy(lineage, metadata, 'Trusted')
            allTrustedGenomeIds.update(trustedGenomeIds)
            print '  Trusted genomes in lineage: ' + str(len(trustedGenomeIds))
                 
            # build gene count table
            print '\nBuilding gene count table.'
            countTable = img.countTable(trustedGenomeIds)
            countTable = img.filterTable(trustedGenomeIds, countTable, 0.9*phyloUbiquityThreshold, 0.9*phyloSingleCopyThreshold)
                    
            # identify marker set
            markerGenes = markerset.markerGenes(trustedGenomeIds, countTable, phyloUbiquityThreshold*len(trustedGenomeIds), phyloSingleCopyThreshold*len(trustedGenomeIds))
            print '  Marker genes: ' + str(len(markerGenes))
        
            tigrToRemove = img.identifyRedundantTIGRFAMs(markerGenes)
            markerGenes = markerGenes - tigrToRemove
            print '  Marker genes after filtering redundant TIGRFAMs: ' + str(len(markerGenes))
            
            phyloMarkerGenes[lineage] = markerGenes
        
        # universal marker genes
        universalMarkerGenes = set()
        for markerGenes in phyloMarkerGenes.values():
            universalMarkerGenes.update(markerGenes)
        
        print ''
        print '  Universal marker genes: ' + str(len(universalMarkerGenes))
        
        # get mapping of marker ids to gene ids for each genome
        genesInGenomes = {}
        for genomeId in allTrustedGenomeIds:
            markerIdToGeneIds = {}
            for line in open(img.genomeDir + '/' + genomeId + '/' + genomeId + '.pfam.tab.txt'):
                lineSplit = line.split('\t')
                markerIdToGeneIds[lineSplit[8]] = markerIdToGeneIds.get(lineSplit[8], []) + [lineSplit[0]]
                
            for line in open(img.genomeDir + '/' + genomeId + '/' + genomeId + '.tigrfam.tab.txt'):
                lineSplit = line.split('\t')
                markerIdToGeneIds[lineSplit[6]] = markerIdToGeneIds.get(lineSplit[6], []) + [lineSplit[0]]
                
            genesInGenomes[genomeId] = markerIdToGeneIds
            
        # get HMM for each universal marker gene
        print ''
        print 'Grabbing HMMs for universal marker genes.'
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
                os.system('hmmfetch -o ./data/hmms/' + markerId + '.hmm /srv/whitlam/bio/db/pfam/27/Pfam-A.hmm ' + markerIdToName[markerId])    
            else:
                os.system('cp /srv/whitlam/bio/db/tigrfam/13.0/' + markerId + '.HMM ./data/hmms/' + markerId + '.hmm')      
        
        # align gene sequences and infer gene trees
        print ''
        print 'Aligning gene sequences.'
        hmmer = HMMERRunner('align')
        for counter, markerId in enumerate(universalMarkerGenes):
            statusStr = '  Finished processing %d of %d (%.2f%%) universal marker genes.' % (counter+1, len(universalMarkerGenes), float(counter+1)*100/len(universalMarkerGenes))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            markerSeqFile = './data/seqs/' + markerId + '.faa'
            fout = open(markerSeqFile, 'w')
            for genomeId in allTrustedGenomeIds:
                seqs = readFasta(img.genomeDir + '/' + genomeId + '/' + genomeId + '.genes.faa')
                
                for geneId in genesInGenomes[genomeId].get(markerId, []):
                    fout.write('>' + genomeId + '_' + geneId + '\n')
                    fout.write(seqs[geneId] + '\n')
            fout.close()
            
            hmmer.align('./data/hmms/' + markerId + '.hmm', markerSeqFile, './data/alignments/' + markerId + '.aln.faa', trim=False, outputFormat='Pfam')
            self.maskAlignment('./data/alignments/' + markerId + '.aln.faa', './data/alignments/' + markerId + '.aln.masked.faa')

        sys.stdout.write('\n')
         
    def maskAlignment(self, inputFile, outputFile):
        """Read HMMER alignment in STOCKHOLM format and output masked alignment in FASTA format."""
        # read STOCKHOLM alignment
        seqs = {}
        for line in open(inputFile):
            line = line.rstrip()
            if line == '' or line[0] == '#' or line == '//':
                if 'GC RF' in line:
                    mask = line.split()[1].rstrip()
                continue
            else:
                lineSplit = line.split()
                seqs[lineSplit[0]] = lineSplit[1].upper().replace('.', '-').rstrip()

        # output masked sequences in FASTA format
        fout = open(outputFile, 'w')
        for seqId, seq in seqs.iteritems():
            fout.write('>' + seqId + '\n')
            
            maskedSeq = ''.join([mask[i] for i in xrange(0, len(seq)) if mask[i] == 'x'])
            fout.write(maskedSeq + '\n')
        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify genes suitable for phylogenetic inference.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-x', '--phylo_ubiquity', help='ubiquity threshold for defining phylogenetic marker genes', type=float, default = 0.90)
    parser.add_argument('-y', '--phylo_single_copy', help='single-copy threshold for defining phylogenetic marker genes', type=float, default = 0.90)
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)
    
    args = parser.parse_args()
    
    phylogeneticInferenceGenes = PhylogeneticInferenceGenes()
    phylogeneticInferenceGenes.run(args.phylo_ubiquity, args.phylo_single_copy, args.threads)
