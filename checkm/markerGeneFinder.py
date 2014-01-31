###############################################################################
#
# markerGeneFinder.py - identify marker genes in genome bins
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
import multiprocessing as mp
import logging

from common import binIdFromFilename, makeSurePathExists

import defaultValues

from hmmer import HMMERRunner, HMMERModelExtractor
from prodigal import ProdigalRunner

class MarkerGeneFinder():
    """Identify marker genes within binned sequences using Prodigal and HMMER."""
    def __init__(self, threads):  
        self.logger = logging.getLogger()                  
        self.totalThreads = threads
        
        self.TAXONOMIC_MARKER_FILE = 1
        self.TREE_MARKER_FILE = 2
        self.HMM_MODELs_FILE = 3

    def find(self, binFiles, outDir, tableOut, hmmerOut, markerFile):
        """Identify marker genes in each bin using prodigal and HMMER."""
        
        # determine type of marker file
        markerFileType = self.__markerFileType(markerFile)
        
        # create HMM model file
        if markerFileType == self.TAXONOMIC_MARKER_FILE:
            hmmModelFile = os.path.join(outDir, 'storage', 'taxonomic.hmm')
            self.__getTaxonomicMarkerHMMs(markerFile, hmmModelFile)
        elif markerFileType == self.TREE_MARKER_FILE:
            pass
        else:
            hmmModelFile = markerFile

        # process each fasta file
        self.threadsPerSearch = max(1, int(self.totalThreads / len(binFiles)))
        self.logger.info("  Identifying marker genes in %d bins with %d threads:" % (len(binFiles), self.totalThreads))
        
        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for binFile in binFiles:
            workerQueue.put(binFile)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__processBin, args = (outDir, tableOut, hmmerOut, hmmModelFile, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportProgress, args = (len(binFiles), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put(None)
        writeProc.join()
        
    def __markerFileType(self, markerFile):
        """Determine type of marker file."""
        with open(markerFile, 'r') as f:
            header = f.readline()
            
        if '# Taxonomic Marker Set' in header:
            return self.TAXONOMIC_MARKER_FILE
        elif '# Tree Marker Sets' in header:
            return self.TREE_MARKER_FILE
        elif 'HMMER3' in header:
            return self.HMM_MODELs_FILE
        else:
            self.logger.error('Unrecognized file type: ' + markerFile)
            sys.exit()
            
    def __getTaxonomicMarkerHMMs(self, markerFile, outputFile):
        """Create HMM file for taxonomic markers."""
        for line in open(markerFile):
            if line[0] == '#':
                continue
            elif 'LINEAGE\t' in line:
                lineage = line.split('\t')[1].rstrip()
            elif 'GENOME\t' in line:
                numGenomes = int(line.split('\t')[1].rstrip())
            elif 'UBIQUITY\t' in line:
                pass
            elif 'SINGLE_COPY\t' in line:    
                pass
            elif 'COLOCATED_DISTANCE\t' in line:   
                pass
            elif 'COLOCATED_GENOME_PERCENTAGE\t' in line: 
                pass
            else:
                markerSets = eval(line.rstrip())
              
        # get list of marker genes  
        markersGene = []
        for markerSet in markerSets:
            for markerGene in markerSet:
                markersGene.append(markerGene)
            
        # report marker set statistics            
        self.logger.info('  Lineage: %s' % lineage)
        self.logger.info('  # reference genomes: %d' % numGenomes)
        self.logger.info('  # marker genes: %d' % len(markersGene))
        self.logger.info('  # marker sets: %d' % len(markerSets))
        self.logger.info('')
        
        # extract marker genes
        modelExtractor = HMMERModelExtractor(self.totalThreads)
        modelExtractor.extract(defaultValues.HMM_MODELS, markersGene, outputFile)
              
    def __processBin(self, outDir, tableOut, hmmerOut, markerFile, queueIn, queueOut):
        """Thread safe bin processing."""      
        while True:    
            binFile = queueIn.get(block=True, timeout=None) 
            if binFile == None:
                break   
            
            binId = binIdFromFilename(binFile)
            binDir = os.path.join(outDir, binId)
            makeSurePathExists(binDir)
    
            # run Prodigal     
            prodigal = ProdigalRunner(binDir) 
            if not prodigal.areORFsCalled():
                prodigal.run(binFile)
    
            # run HMMER
            hmmer = HMMERRunner()
            tableOutPath = os.path.join(binDir, tableOut)
            hmmerOutPath = os.path.join(binDir, hmmerOut)
            hmmer.search(markerFile, prodigal.aaGeneFile, tableOutPath, hmmerOutPath, '--cpu ' + str(self.threadsPerSearch) + ' --notextw -E 1 --domE 1')
    
            queueOut.put(binId)
            
    def __reportProgress(self, numBins, queueIn):
        """Report number of processed bins."""      
        
        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        
        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins)*100/numBins)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
         
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
            