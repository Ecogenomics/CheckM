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
Test stability of marker set for different named taxonomic groups.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse
import random
import multiprocessing as mp

from checkm.lib.img import IMG
from lib.markerSet import MarkerSet

from numpy import mean, std

class MarkerSetStabilityTest(object):
    def __init__(self):
        self.img = IMG()
        self.markerset = MarkerSet()

    def __processLineage(self, metadata, ubiquityThreshold, singleCopyThreshold, minGenomes, queueIn, queueOut):
        """Assess stability of marker set for a specific named taxonomic group."""
        while True:
            lineage = queueIn.get(block=True, timeout=None) 
            if lineage == None:
                break  
            
            genomeIds = self.img.genomeIdsByTaxonomy(lineage, metadata, 'trusted')
            
            markerGenes = []
            perChange = []
            numGenomesToSelect = int(0.9*len(genomeIds))
            if len(genomeIds) >= minGenomes:  
                # calculate marker set for all genomes in lineage          
                geneCountTable = self.img.geneCountTable(genomeIds)
                markerGenes = self.markerset.markerGenes(genomeIds, geneCountTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
                tigrToRemove = self.img.identifyRedundantTIGRFAMs(markerGenes)

                markerGenes = markerGenes - tigrToRemove

                for _ in xrange(0, 100):
                    # calculate marker set for subset of genomes
                    subsetGenomeIds = random.sample(genomeIds, numGenomesToSelect)
                    geneCountTable = self.img.geneCountTable(subsetGenomeIds)
                    subsetMarkerGenes = self.markerset.markerGenes(subsetGenomeIds, geneCountTable, ubiquityThreshold*numGenomesToSelect, singleCopyThreshold*numGenomesToSelect)
                    tigrToRemove = self.img.identifyRedundantTIGRFAMs(subsetMarkerGenes)
                    subsetMarkerGenes = subsetMarkerGenes - tigrToRemove

                    perChange.append(float(len(markerGenes.symmetric_difference(subsetMarkerGenes)))*100.0 / len(markerGenes))

            if perChange != []:
                queueOut.put((lineage, len(genomeIds), len(markerGenes), numGenomesToSelect, mean(perChange), std(perChange)))
            else:
                queueOut.put((lineage, len(genomeIds), len(markerGenes), numGenomesToSelect, -1, -1))
                
    def __storeResults(self, outputFile, totalLineages, writerQueue):
        """Store results to file."""
        
        fout = open(outputFile, 'w')
        fout.write('Lineage\t# genomes\t# markers\t# sampled genomes\tmean % change\tstd % change\n')

        numProcessedLineages = 0
        while True:
            lineage, numGenomes, numMarkerGenes, numSampledGenomes, meanPerChange, stdPerChange = writerQueue.get(block=True, timeout=None)
            if lineage == None:
                break
                    
            numProcessedLineages += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) lineages.' % (numProcessedLineages, totalLineages, float(numProcessedLineages)*100/totalLineages)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            

            fout.write('%s\t%d\t%d\t%d\t%f\t%f\n' % (lineage, numGenomes, numMarkerGenes, numSampledGenomes, meanPerChange, stdPerChange))

        sys.stdout.write('\n')
            
        fout.close()
        
        
    def run(self, outputFile, ubiquityThreshold, singleCopyThreshold, minGenomes, mostSpecificRank, numThreads):
        """Calculate stability of marker sets for named taxonomic groups."""  
        
        print '  Testing stability of marker sets:'
        
        random.seed(1)
        
        # process each sequence in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        metadata = self.img.genomeMetadata()
        lineages = self.img.lineagesByCriteria(metadata, minGenomes, mostSpecificRank)

        for lineage in lineages:
            workerQueue.put(lineage)

        for _ in range(numThreads):
            workerQueue.put(None)
 
        calcProc = [mp.Process(target = self.__processLineage, args = (metadata, ubiquityThreshold, singleCopyThreshold, minGenomes, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target = self.__storeResults, args = (outputFile, len(lineages), writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None, None, None, None, None))
        writeProc.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="assess stability of marker sets",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('output_file', help='Output file.')
    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 10)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 6)
    parser.add_argument('-t', '--threads', help='Threads to use', type=int, default = 32)

    args = parser.parse_args()

    markerSetStabilityTest = MarkerSetStabilityTest()
    markerSetStabilityTest.run(args.output_file, args.ubiquity, args.single_copy, args.min_genomes, args.most_specific_rank, args.threads)
