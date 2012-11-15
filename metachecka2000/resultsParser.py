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
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
import os 
import argparse
from re import compile as re_compile, split as re_split
import uuid

# MetaChecka2000 imports
import defaultValues 

# other local imports
from simplehmmer.simplehmmer import HMMERRunner, HMMERParser, makeOutputFNs

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Mc2kHmmerResultsParser():
    """This class does the job of parsing through the txt output from a hmmer run"""
    def __init__(self, prefix=''):
        # make the output file names
        (self.txtOut, self.hmmOut) = makeOutputFNs(prefix)
        self.results = {}
        self.qLengths = {}
        self.numQs = 0
    
    def analyseResults(self,
                       directory,
                       hmmFile,
                       eCO=defaultValues.__MC2K_DEFAULT_E_VAL__,
                       lengthCO=defaultValues.__MC2K_DEFAULT_LENGTH__,
                       verbose=False,
                       outFile=''
                       ):
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
            storage = HitManager(folder, eCO, lengthCO)
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

class HitManager():
    """Store all the results for a single bin"""
    def __init__(self, name, eCO, lengthCO):
        self.name = name
        self.markers = {}
        self.eCO = eCO
        self.lengthCO = lengthCO
    
    def vetHit(self, hit):
        """See if this hit meets our exacting standards"""
        # we should first check to see if this hit is spurious
        # evalue is the easiest method
        if hit.full_e_value > self.eCO:
            #print "BAD1:", hit.query_name
            return False
        
        # also from the manual, if the bias is significantly large
        # then we shouldn't trust the hit
        if hit.full_bias > hit.full_score*0.5:
            #print "BAD2:", hit.query_name
            return False
        
        # now we can see if we have a long enough match
        alignment_length = float(hit.ali_to - hit.ali_from)
        length_perc = alignment_length/float(hit.target_length) 
        if length_perc < self.lengthCO:
            #print "BAD3:", hit.query_name, length_perc, self.lengthCO
            return False
        
        # this is presumably a good hit.
        return True
    
    def addHit(self, hit):
        """process this hit and add to the markers if OK
        
        hit is an instance of simplehmmer.HmmerHit
        """
        if self.vetHit(hit):
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

class HMMAligner:
    def __init__(self, prefix=''):
        # make output file names
        (self.txtOut, self.hmmOut) = makeOutputFNs(prefix=prefix)
        (txtOut, self.hmmAlign) = makeOutputFNs(prefix=prefix, mode='align')
        

    def makeAlignments(self,
                       directory,
                       hmm,
                       eCO=defaultValues.__MC2K_DEFAULT_E_VAL__,
                       lengthCO=defaultValues.__MC2K_DEFAULT_LENGTH__,
                       prefix='',
                       verbose=False
                       ):
        """main wrapper for making hmm alignments"""
        # first parse through the hmm txt and work out all the 
        # hits that were made
        # we expect directory to contain a collection of folders
        # names after the original bins
        # we need to use this guy to vet the hits
        HM = HitManager('', eCO, lengthCO)
        hit_lookup = {}
        unique_hits = {}
        for folder in os.listdir(directory):
            # we can now build the hmmer_file_name
            hmmer_file_name = os.path.join(directory, folder, self.txtOut)
            hit_lookup[folder] = {}
            with open(hmmer_file_name, 'r') as hmmer_handle:
                try:
                    HP = HMMERParser(hmmer_handle)
                except:
                    print "Error opening HMM file:", fileName
                    raise
                
                while True:
                    hit = HP.next()
                    if hit is None:
                        break
                    if HM.vetHit(hit): # we only want good hits
                        try:
                            hit_lookup[folder][hit.query_accession].append(hit.target_name)
                        except KeyError:
                            hit_lookup[folder][hit.query_accession] = [hit.target_name]
                        try:
                            unique_hits[hit.query_accession] += 1
                        except KeyError:
                            unique_hits[hit.query_accession] = 1
        
        # make temporary files and write the extracted hmms
        single_hmms = {}
        single_hmm_consensi = {}
        if True:
            HA = HMMERRunner(mode='align', prefix=prefix)
            HF = HMMERRunner(mode='fetch', prefix=prefix)
            self.makeAlignmentModels(HF, hmm, unique_hits.keys(), single_hmms, single_hmm_consensi)
        
        # now we need to go through each of the bins and extract the contigs which we need
        for folder in hit_lookup:
            genes_file_name = os.path.join(directory, folder, 'prodigal_out.faa') # the seqs to align
            final_alignment_file = os.path.join(directory, folder, self.hmmAlign) # where to write the whole shebang to
            self.alignBin(HA, hit_lookup[folder], single_hmms, single_hmm_consensi, genes_file_name, final_alignment_file)
        
        # remove the tmp files
        for file_name in single_hmms:
            os.remove(single_hmms[file_name])
        for file_name in single_hmm_consensi:
            os.remove(single_hmm_consensi[file_name])
            
    def alignBin(self, HA, hits, hmms, consensi, fasta, finalFile):
        """Make alignments for this bin"""
        bin_genes = {}
        with open(fasta, 'r') as fh:
            for cid,seq,qual in self.readfq(fh):
                bin_genes[cid] = seq

        align_tmp_file = os.path.join('/tmp', str(uuid.uuid4()))
        align_seq_file = os.path.join('/tmp', str(uuid.uuid4()))
        os.system('echo %s > %s' % (fasta, align_tmp_file))
        for hit in hits:
            # make a multiple fasta consisting of the nodes here and 
            # the hmm consensus, we overwrite on the first go
            os.system('cat %s > %s' % (consensi[hit], align_seq_file))
            for sequence in hits[hit]:
                os.system("echo '>%s' >> %s" % (sequence, align_seq_file))
                os.system('echo %s >> %s' % (bin_genes[sequence], align_seq_file))
                
            # now we can do the alignment
            # we don't want to overwrite the file
            HA.align(hmms[hit], align_seq_file, align_tmp_file, writeMode='>>')
            
            # end of  this hit
            os.system('echo >> %s' % (align_tmp_file))
            os.system('echo "################################################################################" >> %s' % (align_tmp_file))
            os.system('echo >> %s' % (align_tmp_file))

        os.system('cat %s > %s' % (align_tmp_file, finalFile))
        os.remove(align_seq_file)
        os.remove(align_tmp_file)
        
    def makeAlignmentModels(self, HF, hmm, keyz, hmms, consensi):
        """Make a bunch of temporary hmm files which we'll use when aligning"""
        for key in keyz:
            fetch_filename = os.path.join('/tmp', str(uuid.uuid4()))
            emit_filename = os.path.join('/tmp', str(uuid.uuid4()))
            hmms[key] = fetch_filename
            consensi[key] = emit_filename
            HF.fetch(hmm, key, fetch_filename, emit_filename)
            
    def readfq(self, fp): # this is a generator function
        """https://github.com/lh3"""
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break
                
###############################################################################
###############################################################################
###############################################################################
###############################################################################
