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

__prog_name__ = 'findCircularSeqs.py'
__prog_desc__ = 'identify seqs with near identical start and end motifs'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import difflib

from checkm.lib.seqUtils import readFasta

class FindCircularSeqs(object):
    def __init__(self):
      pass

    def hamming(self, str1, str2):
      diffs = 0
      for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
          diffs += 1
      return diffs

    def matchMotif(self, motif, seq, overlap, motifErrors):
      # check for exact match
      posEndMotif = seq.rfind(motif, max(0, len(seq) - overlap))
      if posEndMotif != -1 or motifErrors == 0:
        return posEndMotif

      # find first match with an acceptable number of mismatchers
      for i in xrange(len(seq)-overlap, len(seq)-len(motif)):
        diff = self.hamming(motif, seq[i:i+len(motif)])
        if diff <= motifErrors:
          return i

      return -1

    def run(self, seqFile, length, trim, motifErrors, overlap, acceptDiff, minLen):
      seqs = readFasta(seqFile)

      for seqId, seq in seqs.iteritems():
        if len(seq) < minLen:
          continue

        trimmedSeq = seq[trim:len(seq)-trim]

        startMotif = trimmedSeq[0:length]

        posEndMotif = self.matchMotif(startMotif, trimmedSeq, overlap, motifErrors)
        if posEndMotif != -1:
          diff = self.hamming(trimmedSeq[0:len(trimmedSeq)-posEndMotif], trimmedSeq[posEndMotif:])
          if diff <= acceptDiff:
            print '[Putative Circular Sequence]'
            print seqId
            print 'M: ' + startMotif
            print 'S: ' + trimmedSeq[0:len(trimmedSeq)-posEndMotif]
            print 'E: ' + trimmedSeq[posEndMotif:]
            print 'Motif length: %d' % (len(trimmedSeq)-posEndMotif)
            print 'Hamming distance: % d' % diff
            print 'Sequence length: %d' % len(trimmedSeq)
            print '--------------'

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('seq_file', help='sequences to search for circularity')
    parser.add_argument('-l', '--length', default=50, help='length of motif to use in search')
    parser.add_argument('-t', '--trim', default=15, help='number of bases to trim from start and end of sequence')
    parser.add_argument('-e', '--errors', default=1, help='acceptable mismatched to motif')
    parser.add_argument('-o', '--overlap', default=300, help='search range of start motif at end of sequence')
    parser.add_argument('-d', '--diff', default=2, help='acceptable number of difference between start and end')
    parser.add_argument('-m', '--min_len', default=1000, help='minimum sequences length')

    args = parser.parse_args()

    try:
        findCircularSeqs = FindCircularSeqs()
        findCircularSeqs.run(args.seq_file, args.length, args.trim, args.errors, args.overlap, args.diff, args.min_len)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
