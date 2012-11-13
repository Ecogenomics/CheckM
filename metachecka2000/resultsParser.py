#!/usr/bin/env python
###############################################################################
#                                                                             #
#    resultsParser.py                                                         #
#                                                                             #
#    Parse the output from a run of 'build'                                   #
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
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
import os 
import argparse
from re import compile as re_compile, split as re_split

# MetaChecka2000 imports

# other local imports
from simplehmmer.simplehmmer import HMMERParser, makeOutputFNs

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class MC2KHmmerResultsParser():
    """This class does the job of parsing through the txt output from a hmmer run"""
    def __init__(self, prefix=''):
        # make the output file names
        (self.txtOut, self.hmmOut) = makeOutputFNs(prefix)
        self.results = {}
        self.qLengths = {}
        self.numQs = 0
    
    def analyseResults(self, directory, hmmFile, eCO=1e-10, lengthCO=0.7, verbose=False, outFile=''):
        """Parse the results in the output directory"""
        old_stdout = sys.stdout
        if("" != outFile):
            try:
                # redirect stdout to a file
                sys.stdout = open(outFile, 'w')
            except:
                print "Error diverting stout to file:", outFile, exc_info()[0]
                raise

        # parse the hmm file itself so we can determine the length
        # of the queries
        name_re = re_compile('^NAME')
        length_re = re_compile('^LENG')
        name = ''
        with open(hmmFile, 'r') as hmm_fh:
            for line in hmm_fh:
                if name == '':
                    if name_re.match(line):
                        name = re_split( r'\s+', line.rstrip() )[1] 
                else:
                    if length_re.match(line):
                        self.qLengths[name] = re_split( r'\s+', line.rstrip() )[1]
                        name = ''
                    
        self.numQs = len(self.qLengths)

        # we expect directory to contain a collection of folders
        # names after the original bins
        for folder in os.listdir(directory):
            # somewhere to store results
            storage = MC2KHitStorage(folder, eCO, lengthCO)
            # we can now build the hmmer_file_name
            hmmer_file_name = os.path.join(directory, folder, self.txtOut)
            # and then we can parse it
            self.parseHmmerResults(hmmer_file_name, storage, verbose)
            self.results[folder] = storage
        
        if not verbose:
            self.printHeader()
        for fasta in self.results:
            self.results[fasta].printSummary(self.qLengths, verbose=verbose)

        # restore stdout        
        if("" != outFile):
            try:
                # redirect stdout to a file
                sys.stdout = old_stdout
            except:
                print "Error restoring stdout", exc_info()[0]
                raise
        
        
    def parseHmmerResults(self, fileName, storage, verbose):
        """Parse through a hmmer file and see what's what"""
        with open(fileName, 'r') as hmmer_handle:
            try:
                HP = HMMERParser(hmmer_handle)
            except:
                print "Error opening HMM file:", fileName
                raise
            
            while True:
                hit = HP.next()
                if hit is None:
                    break
                storage.addHit(hit)

    def printHeader(self):
        """Print the NON_VERBOSE header"""
        # keep count of single, double, triple genes etc...
        print "\t".join(['Bin_name','0','1','2','3','4','5+','comp','cont'])

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class MC2KHitStorage():
    """Store all the results for a single bin"""
    def __init__(self, name, eCO, lengthCO):
        self.name = name
        self.markers = {}
        self.eCO = eCO
        self.lengthCO = lengthCO
    
    def addHit(self, hit):
        """process this hit and add to the markers if OK
        
        hit is an instance of simplehmmer.HmmerHit
        """
        # we should first check to see if this hit is spurious
        # evalue is the easiest method
        if hit.full_e_value > self.eCO:
            #print "BAD1:", hit.query_name
            return
        
        # also from the manual, if the bias is significantly large
        # then we shouldn't trust the hit
        if hit.full_bias > hit.full_score*0.5:
            #print "BAD2:", hit.query_name
            return
        
        # now we can see if we have a long enough match
        alignment_length = float(hit.ali_to - hit.ali_from)
        length_perc = alignment_length/float(hit.target_length) 
        if length_perc < self.lengthCO:
            #print "BAD3:", hit.query_name, length_perc, self.lengthCO
            return
        
        # this is presumably a good hit.
        #print hit.query_name, hit.full_score, length_perc, hit.full_e_value
        try:
            self.markers[hit.query_name] += 1
        except KeyError:
            self.markers[hit.query_name] = 1

    def printSummary(self, queries, verbose=False):
        """print out some information about this bin"""
        # if verbose, then just dump the whole thing
        if verbose:
            print "--------------------"
            print self.name
            for marker in queries:
                try:
                    print "%s\t%d" % (marker, self.markers[marker])
                except KeyError:
                    print "%s\t0" % (marker)
                    
            print "TOTAL:\t%d / %d (%0.2f" % (len(self.markers),
                                              len(queries),
                                              100*float(len(self.markers))/float(len(queries))
                                              )+"%)"
            return

        gene_counts = [0]*6
        for marker in queries:
            # we need to limit it form 0 to 5+
            try:
                if self.markers[marker] > 5:
                    marker_count = 5
                else:
                    marker_count = self.markers[marker]
            except KeyError:
                marker_count = 0
            
            gene_counts[marker_count] += 1
        perc_comp = 100*float(len(self.markers))/float(len(queries))
        perc_cont = 100*float(gene_counts[2] + gene_counts[3] + gene_counts[4] + gene_counts[5])/float(len(queries))  
        print "%s\t%s\t%0.2f\t%0.2f" % (self.name,
                                        "\t".join([str(gene_counts[i]) for i in range(6)]),
                                        perc_comp,
                                        perc_cont
                                        )
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################
