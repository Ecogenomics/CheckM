###############################################################################
#
# hmmer.py - runs HMMER and provides functions for parsing output  
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

import os
import sys
import subprocess

import numpy as np

import defaultValues

from common import makeSurePathExists, checkFileExists

class ProdigalError(BaseException): pass

class ProdigalRunner():
    """Wrapper for running HMMER3"""
    def __init__(self):
        # make sure prodigal is installed
        self.checkForProdigal()
        
    def run(self, query, outputDir, translationTable=11):
        makeSurePathExists(outputDir)
        
        aaFile = os.path.join(outputDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_AA__)
        ntFile = os.path.join(outputDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_NT__)
        gffFile = os.path.join(outputDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_GFF__)
        
        cmd = ('prodigal -q -c -m -f gff -g %s -a %s -d %s -i %s > %s' % (str(translationTable), aaFile, ntFile, query, gffFile))
        os.system(cmd)
        
        return aaFile

    def checkForProdigal(self):
        """Check to see if Prodigal is on the system before we try to run it.
    
        We assume that a successful prodigal -h returns 0 and anything
        else returns something non-zero
        """
        try:
            subprocess.call(['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            print "[Error] Make sure prodigal is on your system path."
            sys.exit()

class ProdigalFastaParser():
    """Parses prodigal FASTA output."""
    def __init__(self):
        pass
    
    def genePositions(self, filename):
        checkFileExists(filename)
        
        gp = {}
        for line in open(filename):
            if line[0] == '>':
                lineSplit = line[1:].split()

                geneId = lineSplit[0]
                startPos = int(lineSplit[2])
                endPos = int(lineSplit[4])
        
                gp[geneId] = [startPos, endPos]
        
        return gp
    
class ProdigalGeneFeatureParser():
    """Parses prodigal FASTA output."""
    def __init__(self, filename):
        checkFileExists(filename)
        
        self.genes = {}
        self.lastCodingBase = {}
        
        self.__parseGFF(filename)
        
        self.codingBaseMasks = {}   
        for seqId in self.genes:
            self.codingBaseMasks[seqId] = self.__buildCodingBaseMask(seqId)

    def __parseGFF(self, filename):
        """Parse genes from GFF file."""
        for line in open(filename):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId not in self.genes:
                geneCounter = 0
                self.genes[seqId] = {}
                self.lastCodingBase[seqId] = 0
                
            geneId = seqId + '_' + str(geneCounter)
            geneCounter += 1
            
            start = int(lineSplit[3])
            end = int(lineSplit[4])
  
            self.genes[seqId][geneId] = [start, end]
            self.lastCodingBase[seqId] = max(self.lastCodingBase[seqId], end)
            
    def __buildCodingBaseMask(self, seqId):
        """Build mask indicating which bases in a sequences are coding."""
        
        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        codingBaseMask = np.zeros(self.lastCodingBase[seqId])
        for pos in self.genes[seqId].values():
            codingBaseMask[pos[0]:pos[1]+1] = 1    
            
        return codingBaseMask
    
    def codingBases(self, seqId, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end)."""
        
        # check if sequence has any genes
        if seqId not in self.genes:
            return 0
        
        # set end to last coding base if not specified
        if end == None:
            end = self.lastCodingBase[seqId]
                         
        return np.sum(self.codingBaseMasks[seqId][start:end])
        