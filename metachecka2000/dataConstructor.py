#!/usr/bin/env python
###############################################################################
#                                                                             #
#    dataConstructor                                                          #
#                                                                             #
#    Wrap prodical and hmmer to make data!                                    #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#                888b     d888  .d8888b.   .d8888b.  888    d8P               #  
#                8888b   d8888 d88P  Y88b d88P  Y88b 888   d8P                #  
#                88888b.d88888 888    888        888 888  d8P                 #
#                888Y88888P888 888             .d88P 888d88K                  #
#                888 Y888P 888 888         .od888P"  8888888b                 #
#                888  Y8P  888 888    888 d88P"      888  Y88b                #
#                888   "   888 Y88b  d88P 888"       888   Y88b               #
#                888       888  "Y8888P"  888888888  888    Y88b              #
#                                                                             #
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

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort", "Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.2.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
from os import system, listdir, mkdir
import os.path
import threading
import time

# MetaChecka2000 imports
import defaultValues

# other local imports
from simplehmmer.simplehmmer import HMMERRunner, checkForHMMER, makeSurePathExists

from cogent import DNA, PROTEIN
from cogent.core.genetic_code import GeneticCodes
from cogent.parse.fasta import MinimalFastaParser


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ProdigalError(BaseException): pass
class HMMERError(BaseException): pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Mc2kHmmerDataConstructor():
    """This class runs prodigal and hmmer creating data for parsing"""
    def __init__(self, checkHmmer=True, checkProdigal=True, threads=1):
        if checkHmmer:
            checkForHMMER()
        if checkProdigal:
            checkForProdigal()

        # thready stuff            
        self.varLock = threading.Lock() # we don't have many variables here, so use one lock for everything    
        self.totalThreads = threads
        self.threadPool = threading.BoundedSemaphore(self.totalThreads)

        self.num_files_total = 0
        self.num_files_started = 0
        self.num_files_parsed = 0

    def buildData(self, inFiles, outFolder, hmm, closed, prefix, verbose=False):
        """Main wrapper used for building datasets"""

        if not os.path.exists(outFolder):
            mkdir(outFolder)

        # get a listing of all the files in the directory
        self.num_files_total = len(inFiles)
        if verbose:
            print "Will process %d files" % self.num_files_total
        
        for fasta in inFiles:
            t = threading.Thread(target=self.processFasta,
                                 args=(fasta,outFolder,hmm,closed,prefix,verbose)
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
            
              
    def processFasta(self, fasta, outFolder, hmm, closed, prefix, verbose=False):
        """Thread safe fasta processing"""
        self.threadPool.acquire()
        file_num = 0
        try:
            self.varLock.acquire()
            try:
                self.num_files_started += 1
                file_num = self.num_files_started
                print "Processing file %s (%d of %d)" % (fasta, self.num_files_started, self.num_files_total) 
            finally:
                self.varLock.release()
            
            # make sure we have somewhere to write files to
            out_dir = os.path.join(outFolder, os.path.basename(fasta))
            makeSurePathExists(out_dir)
            
            out_file = open(os.path.join(out_dir,
                                 defaultValues.__MC2K_DEFAULT_TRANSLATE_FILE__),
                                 'w')

            # translate the fasta file in all six frames
            contig_file = open(fasta)
            for name, seq in MinimalFastaParser(contig_file):
                dna_seq = DNA.makeSequence(seq)
                dna_seq.Name = name
                translations = GeneticCodes[11].sixframes(dna_seq)
                frame_number = 0
                for frame in translations:
                    p = PROTEIN.makeSequence(frame, name+'_'+str(frame_count+=1))
                    out_file.write(p.toFasta())

            out_file.close()
            # run prodigal
            #prod_fasta = self.runProdigal(fasta,
            #                              out_dir,
            #                              closed=closed,
            #                              verbose=verbose
            #                              )
            
            
            # run HMMER
            if verbose:
                self.varLock.acquire()
                try:
                    print "    Running hmmer on file %s" % fasta
                finally:
                    self.varLock.release()
                
            HR = HMMERRunner(prefix=prefix)
            HR.search(hmm, out_file.name, out_dir)

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
    
    def runProdigal(self, fasta, outFolder, verbose=False, closed=False):
        """Wrapper for running prodigal"""
        if verbose:
            self.varLock.acquire()
            try:
                print "    Running prodigal on file %s" % fasta
            finally:
                self.varLock.release()
        # make file names
        
        prod_file = os.path.join(outFolder,
                                 defaultValues.__MC2K_DEFAULT_PROD_FN__)
        # work out if we're closing the ends
        if closed:
            cs = "-c"
        else:
            cs = ""
        prod_result = system("prodigal -a %s -i %s -q %s > /dev/null" % (prod_file, fasta, cs))

        if prod_result == 2560:
            # need to rerun with the -p option!
            print "rerunning prodigal with the -p option for %s" % fasta
            prod_result = system("prodigal -a %s -i %s -p meta -q %s > /dev/null" % (prod_file, fasta, cs))

        if prod_result != 0:
            raise ProdigalError("Error running prodigal %d" % prod_result)
        
        return prod_file
    
def checkForProdigal():
    """Check to see if Prodigal is on the system before we try fancy things
    
    We assume that a successful prodigal -h returns 0 and anything
    else returns something non-zero
    """
    # redirect stdout so we don't get mess!
    try:
        exit_status = system('prodigal -h 1> /dev/null 2> /dev/null')
    except:
      print "Unexpected error!", sys.exc_info()[0]
      raise
  
    if exit_status != 0:
        raise ProdigalError("Error attempting to run prodigal, is it in your path?")
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
