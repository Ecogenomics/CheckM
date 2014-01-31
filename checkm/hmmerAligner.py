###############################################################################
#
# hmmerAlign.py - runs HMMER and provides functions for parsing output  
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
import uuid
import logging
import multiprocessing as mp
from collections import defaultdict

import defaultValues

from common import makeSurePathExists
from seqUtils import readFasta
from hmmer import HMMERRunner
from resultsParser import ResultsParser

class HmmerAligner:
    def __init__(self, threads):   
        self.logger = logging.getLogger() 
        self.totalThreads = threads
           
        self.outputFormat = 'Pfam'
        
    def makeAlignments(self,
                       outDir,
                       hmmTableFile,
                       hmmModelFile,
                       bIgnoreThresholds,
                       evalueThreshold,
                       lengthThreshold,
                       alignOutputDir,
                       bSkipOrfCorrection = False,
                       bBestHitOnly = False,
                       bMultipleHitsOnly = False
                       ):
        """Align HMM hits to model."""
        
        self.logger.info("  Extracting marker genes to align.")

        resultsParser = ResultsParser()
        
        # parse HMM information
        resultsParser.parseHmmerModels(hmmModelFile)
        
        # get HMM hits to each bin 
        resultsParser.parseBinHits(outDir, hmmTableFile, bSkipOrfCorrection, bIgnoreThresholds, evalueThreshold, lengthThreshold)
    
        # extract the ORFs to align
        markerSeqs = self.__extractMarkerSeqs(outDir, resultsParser, bBestHitOnly, bMultipleHitsOnly)
        
        # generate individual HMMs required to create multiple sequence alignments
        hmmModelFiles = {}            
        self.__makeAlignmentModels(hmmModelFile, resultsParser.models, hmmModelFiles)
         
        # align each of the marker genes     
        makeSurePathExists(alignOutputDir)
        self.__alignMarkerGenes(markerSeqs, hmmModelFiles, alignOutputDir)
               
        # remove the temporary HMM model files
        for fileName in hmmModelFiles:
            os.remove(hmmModelFiles[fileName])
            
    def __alignMarkerGenes(self, markerSeqs, hmmModelFiles, alignOutputDir):
        """Align marker genes with HMMs in parallel."""
        
        self.logger.info("  Aligning %d marker genes with %d threads:" % (len(markerSeqs), self.totalThreads))

        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for markerId in markerSeqs:
            workerQueue.put(markerId)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        calcProc = [mp.Process(target = self.__alignMarker, args = (markerSeqs, alignOutputDir, hmmModelFiles, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportAlignmentProgress, args = (len(markerSeqs), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put(None)
        writeProc.join()

    def __alignMarker(self, markerSeqs, alignOutputDir, hmmModelFiles, queueIn, queueOut):
        while True:    
            markerId = queueIn.get(block=True, timeout=None) 
            if markerId == None:
                break  
                  
            alignSeqFile = os.path.join(alignOutputDir, markerId + '.unaligned.faa')
            fout = open(alignSeqFile, 'w')
            numSeqs = 0
            for binId, seqs in markerSeqs[markerId].iteritems():
                for seqId, seq in seqs.iteritems():                    
                    fout.write('>' + binId + defaultValues.SEQ_CONCAT_CHAR + seqId + '\n')
                    fout.write(seq + '\n')
                    numSeqs += 1
            fout.close()
            
            if numSeqs > 0:
                outFile = os.path.join(alignOutputDir, markerId + '.aligned.faa')
                HA = HMMERRunner(mode='align')  
                HA.align(hmmModelFiles[markerId], alignSeqFile, outFile, writeMode='>', outputFormat=self.outputFormat, trim=False)  
                
                makedSeqFile = os.path.join(alignOutputDir, markerId + '.masked.faa')
                self.__maskAlignment(outFile, makedSeqFile)
                
            queueOut.put(markerId)
            
    def __reportAlignmentProgress(self, numMarkers, queueIn):
        """Report number of processed bins."""      
        
        numProcessedGenes = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished aligning %d of %d (%.2f%%) marker genes.' % (numProcessedGenes, numMarkers, float(numProcessedGenes)*100/numMarkers)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        
        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedGenes += 1
                statusStr = '    Finished aligning %d of %d (%.2f%%) marker genes.' % (numProcessedGenes, numMarkers, float(numProcessedGenes)*100/numMarkers)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
         
        if self.logger.getEffectiveLevel() <= logging.INFO:
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
            
    def __extractMarkerSeqs(self, outDir, resultsParser, bBestHitOnly, bMultipleHitsOnly):
        """Extract marker sequences from bin."""
        
        markerSeqs = defaultdict(dict)
        for folder in resultsParser.results:
            # read ORFs for bin 
            aaGeneFile = os.path.join(outDir, folder, defaultValues.PRODIGAL_AA)
            binORFs = readFasta(aaGeneFile)
        
            # extract ORFs hitting a marker
            for markerId, hits in resultsParser.results[folder].markerHits.iteritems():
                markerSeqs[markerId][folder] = {}
                
                hits.sort(key = lambda x: x.full_e_value)
                if bBestHitOnly:  
                    topHit = hits[0]       
                    markerSeqs[markerId][folder][topHit.target_name] = self.__extractSeq(topHit.target_name, binORFs)         
                else:
                    if bMultipleHitsOnly and len(hits) < 2:
                        continue
                    
                    for hit in hits:
                        markerSeqs[markerId][folder][hit.target_name] = self.__extractSeq(hit.target_name, binORFs)
                        
        return markerSeqs
                         
    def __extractSeq(self, seqId, seqs):
        if defaultValues.SEQ_CONCAT_CHAR in seqId:
            seqIds = seqId.split(defaultValues.SEQ_CONCAT_CHAR)
            
            seq = ''
            for seqId in seqIds:
                seq += seqs[seqId]
                
            rtnSeq = seq
        else:
            rtnSeq = seqs[seqId]
            
        if rtnSeq[-1] == '*':
            rtnSeq = rtnSeq[0:-1] # remove final '*' inserted by prodigal
            
        return rtnSeq
             
    def __makeAlignmentModels(self, hmmModelFile, modelKeys, hmmModelFiles):
        """Make temporary HMM files used to create HMM alignments.""" 
        
        self.logger.info("  Extracting %d HMM models with %d threads:" % (len(modelKeys), self.totalThreads))
        
        # process each marker in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for modelId in modelKeys:
            fetchFilename = os.path.join('/tmp', str(uuid.uuid4()))
            hmmModelFiles[modelId] = fetchFilename
            workerQueue.put((modelId, fetchFilename))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        calcProc = [mp.Process(target = self.__extractModel, args = (hmmModelFile, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__reportModelExtractionProgress, args = (len(modelKeys), writerQueue))

        writeProc.start()
        
        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()
            
        writerQueue.put(None)
        writeProc.join()
            
    def __extractModel(self, hmmModelFile, queueIn, queueOut):
        """Extract HMM model."""
        HF = HMMERRunner(mode='fetch')
        
        while True:    
            modelId, fetchFilename = queueIn.get(block=True, timeout=None) 
            if modelId == None:
                break  
        
            HF.fetch(hmmModelFile, modelId, fetchFilename)
            
            queueOut.put(modelId)
            
    def __reportModelExtractionProgress(self, numMarkers, queueIn):
        """Report number of extracted HMMs."""      
        
        numModelsExtracted = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished extracting %d of %d (%.2f%%) HMM models.' % (numModelsExtracted, numMarkers, float(numModelsExtracted)*100/numMarkers)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        
        while True:
            modelId = queueIn.get(block=True, timeout=None)
            if modelId == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numModelsExtracted += 1
                statusStr = '    Finished extracting %d of %d (%.2f%%) HMM models.' % (numModelsExtracted, numMarkers, float(numModelsExtracted)*100/numMarkers)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
         
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
    