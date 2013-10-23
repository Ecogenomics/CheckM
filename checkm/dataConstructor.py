#!/usr/bin/env python

###############################################################################
#
# dataConstructor.py - 
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

import sys
from os import mkdir
import os.path
import threading
import time

import defaultValues

from simplehmmer.simplehmmer import HMMERRunner, checkForHMMER, makeSurePathExists

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class HMMERError(BaseException): 
    pass


class HmmerDataConstructor():
    """This class runs prodigal and hmmer creating data for parsing"""
    def __init__(self, checkHmmer=True, threads=1):
        if checkHmmer:
            checkForHMMER()

        # thready stuff            
        self.varLock = threading.Lock() # we don't have many variables here, so use one lock for everything    
        self.totalThreads = threads
        self.threadPool = threading.BoundedSemaphore(self.totalThreads)

        self.num_files_total = 0
        self.num_files_started = 0
        self.num_files_parsed = 0


    def buildData(self, inFiles, outFolder, hmm, prefix, verbose=False, quiet=False):
        """Main wrapper used for building datasets"""

        if not os.path.exists(outFolder):
            mkdir(outFolder)
        # get a listing of all the files in the directory
        self.num_files_total = len(inFiles)
        if verbose:
            print "Will process %d files" % self.num_files_total
        
        for fasta in inFiles:
            t = threading.Thread(target=self.processFasta,
                                 args=(fasta,outFolder,hmm,prefix,verbose, quiet)
                                 )
            #t.daemon = True
            t.start()
            
        while True:
            # don't exit till we're done!
            self.varLock.acquire()
            doned = False
            try:
                doned = self.num_files_parsed >= self.num_files_total
            finally:
                self.varLock.release()
            
            if doned:
                break
            else:
                time.sleep(1)
            
    def translate_six_frames(self, bioseq, table=11):
        revseq = bioseq.reverse_complement()
        for i in range(3):
            yield bioseq[i:].translate(table)
            yield revseq[i:].translate(table)
              
    def processFasta(self, fasta, outFolder, hmm, prefix, verbose=False, quiet=False):
        """Thread safe fasta processing"""
        self.threadPool.acquire()
        try:
            self.varLock.acquire()
            try:
                self.num_files_started += 1
                if not quiet:
                    print "Processing file %s (%d of %d)" % (fasta, self.num_files_started, self.num_files_total) 
            finally:
                self.varLock.release()
            
            # make sure we have somewhere to write files to
            out_dir = os.path.join(outFolder, os.path.basename(fasta))
            makeSurePathExists(out_dir)
            try:
                out_file = open(os.path.join(out_dir,
                                 defaultValues.__CHECKM_DEFAULT_TRANSLATE_FILE__),
                                 'w')
            except:
                print "Cannot open file:", os.path.join(out_dir,
                                 defaultValues.__CHECKM_DEFAULT_TRANSLATE_FILE__)
                sys.exit(1)

            # translate the fasta file in all six frames
            contig_file = open(fasta)
            for rec in SeqIO.parse(contig_file, 'fasta'):
                frame_number = 0
                for frame in self.translate_six_frames(rec.seq):
                    frame_number+=1
                    SeqIO.write(SeqRecord(frame, description='', id=rec.id+"_"+str(frame_number)), out_file, 'fasta')

            out_file.close()

            # run HMMER
            if verbose:
                self.varLock.acquire()
                try:
                    print "    Running hmmer on file %s" % fasta
                finally:
                    self.varLock.release()
                
            HR = HMMERRunner(prefix=prefix)
            HR.search(hmm, out_file.name, out_dir, '--cpu 1')

            # let the world know we've parsed this file
            self.varLock.acquire()
            try:
                self.num_files_parsed += 1
            finally:
                self.varLock.release()
            
        finally:
            # make sure to close the input and output files no matter what
            out_file.close()
            contig_file.close()

            self.threadPool.release()
