###############################################################################
#
# merger.py - identify bins with complementary sets of marker genes
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
import logging

from common import checkDirExists
from resultsParser import ResultsParser

class Merger():
    def __init__(self):
        self.logger = logging.getLogger()
    
    def run(self, binFiles, outDir, hmmTableFile, 
                binIdToHmmModelFile, binIdToMarkerSet, 
                minDeltaComp, maxDeltaCont, 
                minMergedComp, maxMergedCont):   
        checkDirExists(outDir)
        
        self.logger.info('  Comparing marker sets between all pairs of bins.')
        
        # ensure all bins are using the same marker set
        for _, markerSetI in binIdToMarkerSet.iteritems():
            for _, markerSetJ in binIdToMarkerSet.iteritems():
                if markerSetI != markerSetJ:
                    self.logger.error('  [Error] All bins must use the same marker set to assess potential mergers.')
                    sys.exit(0)
            
        # parse HMM information
        resultsParser = ResultsParser()
        resultsParser.parseHmmerModels(binIdToHmmModelFile)
        
        # get HMM hits to each bin 
        resultsParser.parseBinHits(outDir, hmmTableFile)
        
        # determine union and intersection of marker sets for each pair of bins
        outputFile = os.path.join(outDir, "merger.tsv")
        fout = open(outputFile, 'w')
        fout.write('Bin Id 1\tBin Id 2')
        fout.write('\t# markers in bin 1\tBin 1 completeness\tBin 1 contamination')
        fout.write('\t# markers in bin 2\tBin 2 completeness\tBin 2 contamination')
        fout.write('\t# unique markers\t# duplicate markers')
        fout.write('\tDelta completeness\tDelta contamination\tMerger delta')
        fout.write('\tMerged completeness\tMerged contamination\n')
                    
        binMarkerHits = resultsParser.results
        binIds = sorted(binMarkerHits.keys())
        for i in xrange(0, len(binMarkerHits)):
            binIdI = binIds[i]
            markersI = set(binMarkerHits[binIdI].markerHits.keys())
            
            geneCountsI = binMarkerHits[binIdI].geneCounts(binIdToMarkerSet[binIdI], 
                                                           binMarkerHits[binIdI].markerHits, 
                                                           resultsParser.models,
                                                           True)
            completenessI, contaminationI = geneCountsI[6:8]
             
            for j in xrange(i+1, len(binMarkerHits)):               
                binIdJ = binIds[j]
                markersJ = (binMarkerHits[binIdJ].markerHits.keys())
                
                geneCountsJ = binMarkerHits[binIdJ].geneCounts(binIdToMarkerSet[binIdJ], 
                                                               binMarkerHits[binIdJ].markerHits, 
                                                               resultsParser.models,
                                                               True)
                completenessJ, contaminationJ = geneCountsJ[6:8]
                
                # merge together hits from both bins and calculate completeness and contamination
                mergedHits = {}
                for markerId, hits in binMarkerHits[binIdI].markerHits.iteritems():
                    mergedHits[markerId] = list(hits)
                    
                for markerId, hits in binMarkerHits[binIdJ].markerHits.iteritems():
                    if markerId in mergedHits:
                        mergedHits[markerId].extend(hits)
                    else:
                        mergedHits[markerId] = hits
                
                geneCountsMerged = binMarkerHits[binIdI].geneCounts(binIdToMarkerSet[binIdJ], 
                                                                    mergedHits, 
                                                                    resultsParser.models,
                                                                    True)
                completenessMerged, contaminationMerged = geneCountsMerged[6:8]
                
                if not (completenessMerged >= minMergedComp and contaminationMerged < maxMergedCont):
                    continue
                     
                # calculate merged statistics           
                numUnion = len(markersI.union(markersJ))
                numIntersection = len(markersI.intersection(markersJ))
                deltaComp = completenessMerged - max(completenessI, completenessJ)
                deltaCont = contaminationMerged - max(contaminationI, contaminationJ)
                delta = deltaComp - deltaCont 
                
                if deltaComp >= minDeltaComp and deltaCont < maxDeltaCont:    
                    fout.write('%s\t%s\t%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % 
                                                                        (binIdI, binIdJ, 
                                                                         len(markersI), completenessI, contaminationI,
                                                                         len(markersJ), completenessJ, contaminationJ,
                                                                         numUnion, numIntersection, 
                                                                         deltaComp, deltaCont, delta, 
                                                                         completenessMerged, contaminationMerged))
        
        fout.close()
        
        return outputFile
                