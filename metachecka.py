#!/usr/bin/env python
###############################################################################
#
#    metaChecka
#    
#    Given a phylosift directory structure, parse the results and give the 
#    user some informatve information about the bin and eny possible 
#    contaimination that is might contain
#
#    Copyright (C) 2012 Connor Skennerton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

from optparse import OptionParser
import sys
from pprint import pprint
import os
import re
from cogent.parse.fasta import  MinimalFastaParser

###############################################################################
# CODE HERE
###############################################################################
class Marker:
    def __init__(self, marker_name):
        self.name = marker_name
        self.good_count = 0
        self.bad_count = 0

    def __str__(self):
        return "%s\t%s\t%s" % (self.name, self.good_count, self.bad_count)

    def add_good(self):
        self.good_count += 1

    def add_bad(self):
        self.bad_count += 1

    def total_count(self):
        return self.good_count + self.bad_count


class MarkerMapping:
    def __init__(self, contig_name):
        self.contig = contig_name
        self.marker = {}
    
    def add_marker(self, marker_name):
        try:
            self.marker[marker_name] += 1
        except KeyError:
            self.marker[marker_name] = 1


    def __str__(self):
        pass

# reads the output of phylosift and generates counts of 
# the markers and and the contigs they were found in
class MarkerParser(object):
    def __init__(self, directory="./", cutoff=0.15):
        super(MarkerParser, self).__init__()
        self.directory = directory
        self.cutoff = cutoff

    def read_phylosift_directory(self):
        # dict containing the marker name and the contig and count
        marker_contigs = {}
        # dict containing the marker name and the good and bad counts
        marker_counts = {}
        # open up the directory structure and parse the aligned fasta files
        # glob in the fasta files
        for infile in os.listdir(self.directory):
            match = re.search('(PMPROK\d+)\.fasta', infile)
            if match is not None:
                fp = open(self.directory+infile)
                total_seqs = 0
                total_good = 0
                marker_counts[match.group(1)] = Marker(match.group(1))

                for name, seq in MinimalFastaParser(fp):
                    total_seqs+=1
                    num_gaps = 0.0
                    for c in seq:
                        if c == '-':
                            num_gaps+=1.0
                    if num_gaps/float(len(seq)) <= self.cutoff:
                        marker_counts[match.group(1)].add_good()
                        total_good+=1
                    else:
                        marker_counts[match.group(1)].add_bad()
                    
                    try:
                        marker_contigs[name].add_marker(match.group(1))
                    except KeyError, e:
                        marker_contigs[name] = MarkerMapping(name)
                        marker_contigs[name].add_marker(match.group(1))
                #print infile, total_good, total_seqs
        for marker_name, marker_obj in marker_counts.iteritems():
            print marker_obj

        for name, marker in marker_contigs.iteritems():
            print name
            for k, v in marker.marker.iteritems():
                print k, v
        return 0
        

# parse a bam file to generate the length, gc, and coverage 
# information for all of the contigs to plot
class StatsGenerator(object):
    def __init__(self, arg):
        super(StatsGenerator, self).__init__()
        self.arg = arg

# class for holding all of the information about our contigs
# info that will be useful: gc, coverage, length  
class ContigStats(object):
    def __init__(self, name, length, gc, coverage):
        super(ContigStats, self).__init__()
        self.name = name
        self.length = length
        self.gc = gc
        self.coverage = coverage

    def __str__(self):
        return "%s\t%s\t%s\t%s" % (self.name, self.length, self.gc, self.coverage)

    def read_csv(self, file, field_separator="\t"):
        pass
        
def main( options ):
    marker_parser = MarkerParser(options.path)
    marker_parser.read_phylosift_directory()

###############################################################################
# TEMPLATE SUBS
###############################################################################
#
# Entry point, parse command line args and call out to doWork
#
if __name__ == '__main__':
    # intialise the options parser
    parser = OptionParser("\n\n %prog [options]")
    
    # add options here:
    parser.add_option("-i", "--input-directory", dest="path", default='./', 
        help="Input directory containing the output of phylosift align [default: ./]")
    #parser.add_option("-c", "--cat", dest="cat", default=16, help="The cat's age? [default: 16")
    
    # get and check options
    (opts, args) = parser.parse_args()
    
    #
    # compulsory opts
    #
    #if (opts.frog is None ):
    #    print ("**ERROR: %prog : No frog!")
    #    parser.print_help()
    #    sys.exit(1)
    
    # 
    # do what we came here to do
    #
    main(opts)