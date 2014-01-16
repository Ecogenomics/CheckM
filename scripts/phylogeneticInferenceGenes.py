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
from collections import defaultdict
import multiprocessing as mp

from checkm.lib.img import IMG
from lib.markerSet import MarkerSet

from checkm.seqUtils import readFasta
from checkm.hmmer import HMMERRunner

class PhylogeneticInferenceGenes(object):
    def __init__(self):
        pass

    def __getUniversalMarkerGenes(self, phyloUbiquityThreshold, phyloSingleCopyThreshold):
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
        universalMarkerGenes = None
        for markerGenes in phyloMarkerGenes.values():
            if universalMarkerGenes == None:
                universalMarkerGenes = markerGenes
            else:
                universalMarkerGenes.intersection_update(markerGenes)

        fout = open('./data/phylo_marker_set.txt', 'w')
        fout.write(str(universalMarkerGenes))
        fout.close()

        print ''
        print '  Universal marker genes: ' + str(len(universalMarkerGenes))
        
        return allTrustedGenomeIds, universalMarkerGenes
    
    def __genesInGenomes(self, allTrustedGenomeIds):
        genesInGenomes = {}
        for genomeId in allTrustedGenomeIds:
            markerIdToGeneIds = defaultdict(set)
            for line in open(os.path.join(IMG.genomeDir, genomeId, genomeId + IMG.pfamExtension)):
                lineSplit = line.split('\t')
                markerIdToGeneIds[lineSplit[8]].add(lineSplit[0])

            for line in open(os.path.join(IMG.genomeDir, genomeId, genomeId + IMG.tigrExtension)):
                lineSplit = line.split('\t')
                markerIdToGeneIds[lineSplit[6]].add(lineSplit[0])

            genesInGenomes[genomeId] = markerIdToGeneIds
            
        return genesInGenomes
    
    def __fetchMarkerModels(self, universalMarkerGenes):
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
                os.system('hmmfetch -o ./data/phylo_hmms/' + markerId.replace('pfam', 'PF') + '.hmm /srv/whitlam/bio/db/pfam/27/Pfam-A.hmm ' + markerIdToName[markerId] + ' &> /dev/null')
            else:
                os.system('hmmfetch -o ./data/phylo_hmms/' + markerId + '.hmm /srv/whitlam/bio/db/tigrfam/13.0/' + markerId + '.HMM ' + markerId + ' &> /dev/null')
                
    def __alignMarkers(self, allTrustedGenomeIds, universalMarkerGenes, genesInGenomes, numThreads, outputDir):
        """Perform multithreaded alignment of marker genes using HMM align."""
        
        print ''
        print 'Aligning gene sequences.'
        
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        for _, markerId in enumerate(universalMarkerGenes):
            workerQueue.put(markerId)
            
        for _ in range(numThreads):
            workerQueue.put(None)
            
        calcProc = [mp.Process(target = self.__runHmmAlign, args = (allTrustedGenomeIds, genesInGenomes, outputDir, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target = self.__reportThreads, args = (len(universalMarkerGenes), writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put(None)
        writeProc.join()
        
    def __runHmmAlign(self, allTrustedGenomeIds, genesInGenomes, outputDir, queueIn, queueOut):
        """Run each marker gene in a separate thread."""
        
        while True:
            markerId = queueIn.get(block=True, timeout=None) 
            if markerId == None:
                break 
            
            print markerId
            
            modelName = markerId
            if modelName.startswith('pfam'):
                modelName = modelName.replace('pfam', 'PF')

            markerSeqFile = './data/gene_seqs/' + modelName + '.faa'
            fout = open(markerSeqFile, 'w')
            for genomeId in allTrustedGenomeIds:
                seqs = readFasta(IMG.genomeDir + '/' + genomeId + '/' + genomeId + '.genes.faa')

                for geneId in genesInGenomes[genomeId].get(markerId, []):
                    fout.write('>' + genomeId + '_' + geneId + '\n')
                    fout.write(seqs[geneId] + '\n')
            fout.close()
            
            hmmer = HMMERRunner('align')
            hmmer.align(os.path.join('./data/phylo_hmms', modelName + '.hmm'), markerSeqFile, os.path.join(outputDir, modelName + '.aln.faa'), trim=False, outputFormat='Pfam')
            self.__maskAlignment(os.path.join(outputDir, modelName + '.aln.faa'), os.path.join(outputDir, modelName + '.aln.masked.faa'))
            
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

    def run(self, phyloUbiquityThreshold, phyloSingleCopyThreshold, numThreads, outputDir):
        # make sure output directory is empty
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        files = os.listdir(outputDir)
        for f in files:
            os.remove(os.path.join(outputDir, f))

        # get universal marker genes
        allTrustedGenomeIds, universalMarkerGenes = self.__getUniversalMarkerGenes(phyloUbiquityThreshold, phyloSingleCopyThreshold)

        # get mapping of marker ids to gene ids for each genome
        genesInGenomes = self.__genesInGenomes(allTrustedGenomeIds)

        # get HMM for each universal marker gene
        self.__fetchMarkerModels(universalMarkerGenes)

        # align gene sequences and infer gene trees
        self.__alignMarkers(allTrustedGenomeIds, universalMarkerGenes, genesInGenomes, numThreads, outputDir)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify genes suitable for phylogenetic inference.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-x', '--phylo_ubiquity', help='ubiquity threshold for defining phylogenetic marker genes', type=float, default = 0.90)
    parser.add_argument('-y', '--phylo_single_copy', help='single-copy threshold for defining phylogenetic marker genes', type=float, default = 0.90)
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)
    parser.add_argument('-o', '--output_dir', help='output directory for alignments', default = './data/gene_alignments/')

    args = parser.parse_args()

    phylogeneticInferenceGenes = PhylogeneticInferenceGenes()
    phylogeneticInferenceGenes.run(args.phylo_ubiquity, args.phylo_single_copy, args.threads, args.output_dir)
