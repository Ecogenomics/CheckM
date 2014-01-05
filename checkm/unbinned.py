###############################################################################
#
# unbinned.py - identify unbinned sequences
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

from seqUtils import readFasta
import logging

from common import checkFileExists, reassignStdOut, restoreStdOut

class Unbinned():
    def __init__(self):
        self.logger = logging.getLogger()
    
    def run(self, binFiles, seqFile, outFile, minSeqLen):       
        checkFileExists(seqFile)
        
        # get list of sequences in bins
        self.logger.info('  Reading binned sequences.')
            
        binnedSeqs = {}
        for binFile in binFiles:
            seqs = readFasta(binFile)
            binnedSeqs.update(seqs)
            
        self.logger.info('    Read %d binned sequences.' % (len(binnedSeqs)))
            
        # get list of all sequences
        self.logger.info('  Reading all sequences.')
        allSeqs = readFasta(seqFile)
        self.logger.info('    Read %d sequences.' % (len(allSeqs)))
        
        # write all unbinned sequences
        self.logger.info('  Identifying unbinned sequences.')
        
        # redirect output
        oldStdOut = reassignStdOut(outFile)
        
        unbinnedCount = 0
        for seqId, seq in allSeqs.iteritems():
            if seqId not in binnedSeqs:
                if len(seq) >= minSeqLen:
                    unbinnedCount += 1
                    print('>' + seqId)
                    print(seq)
                    
        # restore stdout   
        restoreStdOut(outFile, oldStdOut)
        
        self.logger.info('    Identified %d unbinned sequences.' % (unbinnedCount))
                