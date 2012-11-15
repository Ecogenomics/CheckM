#!/usr/bin/env python
###############################################################################
#
# metaChecka
#
# Given a phylosift directory structure, parse the results and give the
# user some informatve information about the bin and eny possible
# contaimination that is might contain
#
# Copyright (C) 2012 Connor Skennerton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################


import sys
import os
import re
import csv
from optparse import OptionParser
from pprint import pprint
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.usage import BaseUsage
import numpy as np
import matplotlib
matplotlib.use("PS")
import matplotlib.pyplot as plt

###############################################################################
# CODE HERE
###############################################################################
class Marker:
    def __init__(self, marker_name):
        self.name = marker_name
        self.good_count = 0
        self.bad_count = 0
        self.contig = {}

    def __str__(self):
        return "%s\t%s\t%s" % (self.name, self.good_count, self.bad_count)

    def add_good(self, contig_name):
        self.good_count += 1
        self.contig[contig_name] = True

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
                        marker_counts[match.group(1)].add_good(name)
                        total_good+=1
                    else:
                        marker_counts[match.group(1)].add_bad()
                    
                    try:
                        marker_contigs[name].add_marker(match.group(1))
                    except KeyError, e:
                        marker_contigs[name] = MarkerMapping(name)
                        marker_contigs[name].add_marker(match.group(1))
                #print infile, total_good, total_seqs
        #for marker_name, marker_obj in marker_counts.iteritems():
        # print marker_obj

        #for name, marker in marker_contigs.iteritems():
        # print name
        # for k, v in marker.marker.iteritems():
        # print k, v
        return marker_contigs, marker_counts
        
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


# parse a bam file to generate the length, gc, and coverage
# information for all of the contigs to plot
class StatsGenerator(object):
    def __init__(self):
        super(StatsGenerator, self).__init__()
        self.contigs = {}
        #self.arg = arg

    def read_csv(self, csv_file, field_separator=","):
        fd = open(csv_file)
        for line in fd:
            try:
                name, length, gc, coverage = line.rstrip().split(field_separator)
            except ValueError, e:
                print e
                print "Offending line: "
                print line
                raise e
            cs = ContigStats(name, int(length), float(gc), float(coverage))
            self.contigs[name] = cs

    def read_fasta(self, fasta_file, contig_list):
        fd = open(fasta_file)
        for name, seq in MinimalFastaParser(fd):
            base_usage = BaseUsage(seq)
            gc = base_usage.content("GC")
            length = len(seq)

    def read_bam(self, bam_file, contig_list):
        try:
            import pysam
            try:
                samFile = pysam.Samfile(bam_file, 'rb')
            except:
                print "Unable to open BAM file -- did you supply a SAM file instead?"
                sys.exit(1)
            
            basename = os.path.splitext(bam_file)[0]
            output_csv = csv.writer(open(basename + '.cov.csv', "wb"),
                quoting = csv.QUOTE_NONNUMERIC )
            
            for reference, length in zip(samFile.references, samFile.lengths):
                #num_reads = samFile.count(reference, 0, length)
                #print num_reads / length
                tid = samFile.gettid(reference)
                if(tid != -1):
                    tmp_cov = 0
                    cov_vec = []
                    for base in samFile.pileup(samFile.header['SQ'][tid]['SN'],
                                            0, length):
                        if opts.averages:
                            tmp_cov += base.n
                        else:
                            cov_vec.append(base.n)
                    if opts.averages:
                        ave_coverage = float(tmp_cov) / float(length)
                        print reference, ave_coverage
                    else:
                        if len(cov_vec) == 0:
                            cov_vec = [0]*length
                        cov_vec.insert(0,reference)
                        output_csv.writerow(cov_vec)

        except ImportError, e:
            print "ERROR: 'pysam' not in pythonpath, cannot read bam file"
            raise e


# Here we plot the contigs with the markers to see where they all fall
class Plotter(object):
    def __init__(self, contigs, coverage_limit, output="test" ):
            
        super(Plotter, self).__init__()
        self.contigs = contigs
        self.out_name = output
        if coverage_limit is not None:
            self.cov_lim = coverage_limit.split(":")
        else:
            self.cov_lim = None

    def plot_markers(self, markers, marker_mapping):
        gc_list = []
        gc_marker_list = []
        gc_single_marker_list = []
        length_list = []
        length_marker_list = []
        length_single_marker_list = []
        coverage_list = []
        coverage_marker_list = []
        coverage_single_marker_list = []
        
        for k,v in self.contigs.iteritems():
            # check if this contig contains a marker
            if k in marker_mapping:
                for marker_name, marker_count in marker_mapping[k].marker.iteritems():
                    if markers[marker_name].good_count > 1:
                        gc_marker_list.append(v.gc)
                        length_marker_list.append(v.length)
                        coverage_marker_list.append(v.coverage)
                    else:
                        gc_single_marker_list.append(v.gc)
                        length_single_marker_list.append(v.length)
                        coverage_single_marker_list.append(v.coverage)
            else:
                gc_list.append(v.gc)
                length_list.append(v.length)
                coverage_list.append(v.coverage)

        fig = plt.figure()
        
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        if length_list and gc_list:
            ax1.scatter(length_list, gc_list, facecolors='none')
        if length_marker_list and gc_marker_list:
            ax1.scatter(length_marker_list, gc_marker_list, c='r')
        if length_single_marker_list and gc_single_marker_list:
            ax1.scatter(length_single_marker_list, gc_single_marker_list, c='g')
        ax1.set_xlabel('length')
        ax1.set_ylabel('GC %')
        ax1.set_xlim(left=0)
        #ax1.set_ylim(bottom=0)
        
        if length_list and coverage_list:
            ax2.scatter(length_list, coverage_list, facecolors='none')
        if length_marker_list and coverage_marker_list:
            ax2.scatter(length_marker_list, coverage_marker_list, c='r')
        if length_single_marker_list and coverage_single_marker_list:
            ax2.scatter(length_single_marker_list, coverage_single_marker_list,
                    c='g')
        ax2.set_xlabel('length')
        ax2.set_ylabel('coverage')
        ax2.set_xlim(left=0)
        if self.cov_lim is not None:
            try:
                ax2.set_ylim(int(self.cov_lim[0]), int(self.cov_lim[1]))
            except IndexError, e:
                print e
                print "The list is:"
                print self.cov_limit
                raise e
        
        if gc_list and coverage_list:
            ax3.scatter(gc_list, coverage_list, facecolors='none')
        if gc_marker_list and coverage_marker_list:
            ax3.scatter(gc_marker_list, coverage_marker_list, c='r')
        if gc_single_marker_list and coverage_single_marker_list:
            ax3.scatter(gc_single_marker_list, coverage_single_marker_list, c='g')
        ax3.set_xlabel('GC %')
        ax3.set_ylabel('coverage')
        if self.cov_lim is not None:
            ax3.set_ylim(int(self.cov_lim[0]), int(self.cov_lim[1]))

        fig.canvas.draw()
        plt.savefig(self.out_name+".ps")
                        

      
def main( options ):
    # read the phylosift directory to figure out which
    # contigs have which markers on them
    marker_parser = MarkerParser(options.path)
    marker_contigs, marker_counts = marker_parser.read_phylosift_directory()
    
    # read in the contig stats file
    sg = StatsGenerator()
    sg.read_csv(options.contig_stats)
    #for k,v in sg.contigs.iteritems():
    # print v
    plotter = Plotter(sg.contigs, options.cov_lim, options.out_name )
    plotter.plot_markers(marker_counts, marker_contigs)

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
    parser.add_option("-c", "--contig_stats", dest="contig_stats",
        help="A comma separated file containing assembly stats for one contig per line.\
The columns should be: name, length, gc, coverage")
    parser.add_option("-o", "--output", dest="out_name", default="test",
        help="output name prefix for postscript file (no file extension)")
    parser.add_option("-y", "--cov_lim", dest="cov_lim",help="The upper and\
lower limits of the coverage to outputi in the form of 'lower:upper'")
    # get and check options
    (opts, args) = parser.parse_args()
    
    #
    # compulsory opts
    #
    if (opts.path is None or opts.contig_stats is None):
        print ("Please supply a file containing contig stats and the output from the phylosift align directory")
        parser.print_help()
        sys.exit(1)
    
    #
    # do what we came here to do
    #
    main(opts)
