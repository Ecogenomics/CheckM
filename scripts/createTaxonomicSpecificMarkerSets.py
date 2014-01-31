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

__prog_name__ = 'TaxonomicMarkerSets'
__prog_desc__ = 'create all taxonomic-specific marker sets'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import multiprocessing as mp

from checkm.lib.img import IMG
from lib.markerSet import MarkerSet
from lib.taxonomyUtils import rankPrefixes

class TaxonomicMarkerSets(object):
    def __init__(self):
        pass

    def __workerThread(self, ubiquityThreshold, singleCopyThreshold, 
                       minGenomes, 
                       colocatedDistThreshold, colocatedGenomeThreshold, 
                       metadata, 
                       queueIn, queueOut):
        """Process each data item in parallel."""
        
        img = IMG()
        markerset = MarkerSet()

        while True:
            lineage = queueIn.get(block=True, timeout=None)
            if lineage == None:
                break

            genomeIds = img.genomeIdsByTaxonomy(lineage, metadata, 'trusted')
            genomeMissingData = img.genomesWithMissingData(genomeIds)
            genomeIds -= genomeMissingData
            if len(genomeIds) >= minGenomes:
                geneCountTable = img.geneCountTable(genomeIds)
                geneCountTable = img.filterGeneCountTable(genomeIds, geneCountTable, 0.9*ubiquityThreshold, 0.9*singleCopyThreshold)
                
                markerGenes = markerset.markerGenes(genomeIds, geneCountTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
                tigrToRemove = img.identifyRedundantTIGRFAMs(markerGenes)
                markerGenes = markerGenes - tigrToRemove
    
                geneDistTable = img.geneDistTable(genomeIds, markerGenes)
                colocatedGenes = markerset.colocatedGenes(geneDistTable, colocatedDistThreshold, colocatedGenomeThreshold)
                colocatedSets = markerset.colocatedSets(colocatedGenes, markerGenes)
            else:
                colocatedSets = None

            # allow results to be processed or written to file
            queueOut.put((lineage, colocatedSets, len(genomeIds)))

    def __writerThread(self, pfamIdToPfamAcc,
                       ubiquityThreshold, singleCopyThreshold, 
                       colocatedDistThreshold, colocatedGenomeThreshold,
                       outputDir, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            lineage, colocatedSets, numGenomes = writerQueue.get(block=True, timeout=None)
            if lineage == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) lineages.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            if colocatedSets != None:
                taxonomy = lineage.split(';')
                rankPrefix = rankPrefixes[len(taxonomy)-1]
                
                cladeName = taxonomy[-1].strip().replace(' ', '_')
                fout = open(os.path.join(outputDir, rankPrefix + cladeName + '.txt'), 'w')
                fout.write('# Taxonomic Marker Set\n')
                fout.write('LINEAGE\t' + lineage + '\n')
                fout.write('GENOME\t' + str(numGenomes) + '\n')
                fout.write('UBIQUITY\t' + str(ubiquityThreshold) + '\n')
                fout.write('SINGLE_COPY\t' + str(singleCopyThreshold) + '\n')
                fout.write('COLOCATED_DISTANCE\t' + str(colocatedDistThreshold) + '\n')
                fout.write('COLOCATED_GENOME_PERCENTAGE\t' + str(colocatedGenomeThreshold) + '\n')
                
                # change model names to accession numbers, and make
                # sure there is an HMM model for each PFAM
                mungedColocatedSets = []
                for geneSet in colocatedSets:
                    s = set()
                    for geneId in geneSet:
                        if 'pfam' in geneId:
                            pfamId = geneId.replace('pfam', 'PF')
                            if pfamId not in pfamIdToPfamAcc:
                                print geneId # missing PFAM HMM
                            else:                         
                                s.add(pfamIdToPfamAcc[pfamId])
                        else:
                            s.add(geneId)
                    mungedColocatedSets.append(s)
                            
                fout.write(str(mungedColocatedSets))
                fout.close()

        sys.stdout.write('\n')
        
    def __pfamIdToPfamAcc(self, img):
        pfamIdToPfamAcc = {}
        for line in open(img.pfamHMMs):
            if 'ACC' in line:
                acc = line.split()[1].strip()
                pfamId = acc.split('.')[0]
                
                pfamIdToPfamAcc[pfamId] = acc
        
        return pfamIdToPfamAcc
        
    def run(self, outputDir, ubiquityThreshold, singleCopyThreshold, minGenomes, colocatedDistThreshold, colocatedGenomeThreshold, threads):
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
            
        # determine lineages to process
        img = IMG()
        metadata = img.genomeMetadata()
        lineages = img.lineagesSorted(metadata)
        
        # determine HMM model accession numbers
        pfamIdToPfamAcc = self.__pfamIdToPfamAcc(img)
        
        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for lineage in lineages:
            workerQueue.put(lineage)

        for _ in range(threads):
            workerQueue.put(None)

        workerProc = [mp.Process(target = self.__workerThread, args = (ubiquityThreshold, singleCopyThreshold, 
                                                                       minGenomes, 
                                                                       colocatedDistThreshold, colocatedGenomeThreshold, 
                                                                       metadata,
                                                                       workerQueue, writerQueue)) for _ in range(threads)]
        writeProc = mp.Process(target = self.__writerThread, args = (pfamIdToPfamAcc, 
                                                                       ubiquityThreshold, singleCopyThreshold, 
                                                                       colocatedDistThreshold, colocatedGenomeThreshold,
                                                                       outputDir, len(lineages), writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put((None, None, None))
        writeProc.join()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-u', '--ubiquity', help='ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='minimum genomes required to infer marker set', type=int, default = 5)
    parser.add_argument('-d', '--distance_threshold', help='genomic distance to be considered co-located', type=float, default=5000)
    parser.add_argument('-g', '--genome_threshold', help='percentage of genomes required to be considered co-located', type=float, default=0.95)
    parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

    args = parser.parse_args()

    try:
        taxonomicMarkerSets = TaxonomicMarkerSets()
        taxonomicMarkerSets.run(args.output_dir, args.ubiquity, args. single_copy, args.min_genomes, args.distance_threshold, args.genome_threshold, args.threads)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
