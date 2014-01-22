###############################################################################
#
# profile.py - calculate percentage of reads mapped to each bin
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
import logging

import defaultValues
from common import checkFileExists, reassignStdOut, restoreStdOut

class Profile():
    def __init__(self):
        self.logger = logging.getLogger()
    
    def run(self, coverageFile, outFile):
        checkFileExists(coverageFile)
        
        # get number of reads mapped to each bin
        self.logger.info('  Determining number of reads mapped to each bin.')
            
        readsMappedToBin = {}
        binSize = {}
        totalMappedReads = {}
        bHeader = True
        for line in open(coverageFile):
            if bHeader:
                bHeader = False
                continue
            
            lineSplit = line.split('\t')
            
            #seqId = lineSplit[0]
            binId = lineSplit[1]
            
            seqLen = int(lineSplit[2])
            binSize[binId] = binSize.get(binId, 0) + seqLen
            
            if binId not in readsMappedToBin:
                readsMappedToBin[binId] = {}

            for i in xrange(3, len(lineSplit), 3):
                bamId = lineSplit[i]
                mappedReads = int(lineSplit[i+2])        
                
                totalMappedReads[bamId] = totalMappedReads.get(bamId, 0) + mappedReads       
                readsMappedToBin[binId][bamId] = readsMappedToBin[binId].get(bamId, 0) + mappedReads
                
        # calculate percentage of mapped reads and binned populations
        perMappedReads = {}
        normBinCoverage = {}
        sumNormBinCoverage = {}
        for binId, bamIds in readsMappedToBin.iteritems():       
            perMappedReads[binId] = {}
            normBinCoverage[binId] = {}
            
            for bamId in bamIds:
                perMR = float(readsMappedToBin[binId][bamId]) / totalMappedReads[bamId]
                perMappedReads[binId][bamId] = perMR
                
                if binId == defaultValues.UNBINNED:
                    continue
                                
                normCoverage = perMR / binSize[binId]
                normBinCoverage[binId][bamId] = normCoverage 
                sumNormBinCoverage[bamId] = sumNormBinCoverage.get(bamId, 0) + normCoverage
                
        for binId, bamIds in normBinCoverage.iteritems():
            for bamId in bamIds:
                normBinCoverage[binId][bamId] /= sumNormBinCoverage[bamId]

        # write community profile
        oldStdOut = reassignStdOut(outFile)
        
        sortedBinIds = sorted(readsMappedToBin.keys())
        sortedBamIds = sorted(readsMappedToBin[sortedBinIds[0]].keys())
        
        header = 'Bin Id\tBinSize'
        for bamId in sortedBamIds:
            header += '\t' + bamId + ': mapped reads'
            header += '\t' + bamId + ': % mapped reads'
            header += '\t' + bamId + ': % binned populations'
            header += '\t' + bamId + ': % community'
        print(header)
        
        for binId in sortedBinIds:
            rowStr = binId
            rowStr += '\t' + str(binSize[binId])
            
            for bamId in sortedBamIds:
                rowStr += '\t' + str(readsMappedToBin[binId][bamId])
                rowStr += '\t' + str(perMappedReads[binId][bamId]*100.0)
                 
                if defaultValues.UNBINNED in perMappedReads:
                    unbinnedPercentage = perMappedReads[defaultValues.UNBINNED][bamId]
                else:
                    unbinnedPercentage = 0
                    
                if binId == defaultValues.UNBINNED:
                    rowStr += '\t' + 'NA' 
                    rowStr += '\t' + str(unbinnedPercentage*100.0)
                else:
                    rowStr += '\t' + str(normBinCoverage[binId][bamId]*100.0)
                    rowStr += '\t' + str(normBinCoverage[binId][bamId]*100.0 * (1.0 - unbinnedPercentage))
            
            print(rowStr)
        
        restoreStdOut(outFile, oldStdOut)