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
Use 'add-one-more' analysis to assess stability of marker sets for varying 
  numbers of reference genomes.
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

class MarkerSetStability(object):
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
            
            changeMarkerSetSize = {}
            markerGenes = []
            if len(genomeIds) >= minGenomes:  
                # calculate marker set for all genomes in lineage          
                geneCountTable = self.img.geneCountTable(genomeIds)
                markerGenes = self.markerset.markerGenes(genomeIds, geneCountTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
                tigrToRemove = self.img.identifyRedundantTIGRFAMs(markerGenes)
                markerGenes = markerGenes - tigrToRemove
     
                for selectPer in range(50, 101, 5):
                    numGenomesToSelect = int(float(selectPer)/100 * len(genomeIds))
                    perChange = []
                    for _ in range(0, 10):
                        # calculate marker set for subset of genomes
                        subsetGenomeIds = random.sample(genomeIds, numGenomesToSelect)
                        geneCountTable = self.img.geneCountTable(subsetGenomeIds)
                        subsetMarkerGenes = self.markerset.markerGenes(subsetGenomeIds, geneCountTable, ubiquityThreshold*numGenomesToSelect, singleCopyThreshold*numGenomesToSelect)
                        tigrToRemove = self.img.identifyRedundantTIGRFAMs(subsetMarkerGenes)
                        subsetMarkerGenes = subsetMarkerGenes - tigrToRemove
    
                        perChange.append(float(len(markerGenes.symmetric_difference(subsetMarkerGenes)))*100.0 / len(markerGenes))
    
                    changeMarkerSetSize[selectPer] = [mean(perChange), std(perChange)]  

            queueOut.put((lineage, len(genomeIds), len(markerGenes), changeMarkerSetSize))

    def __storeResults(self, outputFile, totalLineages, writerQueue):
        """Store results to file."""
        
        fout = open(outputFile, 'w')
        fout.write('Lineage\t# genomes\t# markers\tsubsample %\tmean % change\tstd % change\n')

        numProcessedLineages = 0
        while True:
            lineage, numGenomes, numMarkerGenes, changeMarkerSetSize = writerQueue.get(block=True, timeout=None)
            if lineage == None:
                break
                    
            numProcessedLineages += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) lineages.' % (numProcessedLineages, totalLineages, float(numProcessedLineages)*100/totalLineages)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            for selectPer in sorted(changeMarkerSetSize.keys()): 
                fout.write('%s\t%d\t%d\t%d\t%f\t%f\n' % (lineage, numGenomes, numMarkerGenes, selectPer, changeMarkerSetSize[selectPer][0], changeMarkerSetSize[selectPer][1]))

        sys.stdout.write('\n')
            
        fout.close()
        
        
    def run(self, outputFile, ubiquityThreshold, singleCopyThreshold, minGenomes, mostSpecificRank, numThreads):
        """Calculate stability of marker sets for named taxonomic groups."""  
        
        print('  Calculating stability of marker sets:')
        
        random.seed(1)
        
        # process each sequence in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        metadata = self.img.genomeMetadata()
        lineages = self.img.lineagesByCriteria(metadata, minGenomes, mostSpecificRank)
        
        #lineages = ['Bacteria']
        #lineages += ['Bacteria;Proteobacteria']
        #lineages += ['Bacteria;Proteobacteria;Gammaproteobacteria']
        #lineages += ['Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales']
        #lineages += ['Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae']
        #lineages += ['Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia']
        #lineages += ['Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;coli']
        
        #lineages = ['Archaea']
        #lineages += ['Archaea;Euryarchaeota']
        #lineages += ['Archaea;Euryarchaeota;Methanomicrobia']
        #lineages += ['Archaea;Euryarchaeota;Methanomicrobia;Methanosarcinales']
        #lineages += ['Archaea;Euryarchaeota;Methanomicrobia;Methanosarcinales;Methanosarcinaceae']

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

        writerQueue.put((None, None, None, None))
        writeProc.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="assess stability of marker sets",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('output_file', help='Output file.')
    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 20)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default = 6)
    parser.add_argument('-t', '--threads', help='Threads to use', type=int, default = 32)

    args = parser.parse_args()

    markerSetStability = MarkerSetStability()
    markerSetStability.run(args.output_file, args.ubiquity, args.single_copy, args.min_genomes, args.most_specific_rank, args.threads)
