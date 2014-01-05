###############################################################################
#
# aminoAcidIdentity.py - calculate AAI between aligned marker genes
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

import multiprocessing as mp
import logging

class AminoAcidIdentity():
    """Calculate coverage of all sequences."""
    def __init__(self, threads):
        self.logger = logging.getLogger()
        
        # thready stuff
        self.totalThreads = threads
        self.threadPool = mp.BoundedSemaphore(self.totalThreads)
        self.manager = mp.Manager()

    def calculate(self, bamFiles, outFile, bPairsOnly=False):
        """Calculate coverage of sequences for each BAM file."""
        pass

    def __calculateResults(self, bamIn, queueIn, queueOut, bPairsOnly):
        """Calculate coverage of reference sequences using multiple processes."""
        pass

    def __storeResults(self, coverageInfo, queue, numReferences):
        """Store coverage information results determined by multiple processes (i.e., to shared memory)."""
        pass

    def __processBam(self, bamFile, bPairsOnly):
        """Calculate coverage of sequences in BAM file."""
        pass

    def __processReference(self, bamIn, reference, length, bPairsOnly):
        """Calculate coverage of reference sequences."""
        pass

