###############################################################################
#
# seqUtils.py - Common functions for interacting with sequences
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

def readFasta(fastaFile):
    seqs = {}
    for line in open(fastaFile):
        if line[0] == '>':
            seqId = line[1:].split()[0]
            seqs[seqId] = ''
        else:
            seqs[seqId] += line.rstrip()
            
    return seqs

def baseCount(seq):
    testSeq = seq.upper()
    a = testSeq.count('A')
    c = testSeq.count('C')
    g = testSeq.count('G') 
    t = testSeq.count('T') + testSeq.count('U')
    
    return a, c, g, t

def calculateN50(seqLens):
    thresholdN50 = sum(seqLens) / 2.0
    
    seqLens.sort(reverse=True)
    
    testSum = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break
        
    return N50