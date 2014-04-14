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

"""
Produce plots comparing performance of different marker sets (best vs. selected vs. domain).
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse

from collections import defaultdict, Counter

from lib.plots.boxplot import BoxPlot
from checkm.lib.img import IMG

from numpy import mean, abs, array

from scipy.stats import scoreatpercentile

class SimComparePlots(object):
    def __init__(self):
        self.img = IMG()
  
    def __readResults(self, filename):
        results = defaultdict(dict)
        genomeIds = set()
        with open(filename) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                simId = lineSplit[0]
                genomeId = simId.split('-')[0]
                genomeIds.add(genomeId)
                
                bestCompIM = [float(x) for x in lineSplit[6].split(',')]
                bestContIM = [float(x) for x in lineSplit[7].split(',')]
                                
                domCompIM = [float(x) for x in lineSplit[8].split(',')]
                domContIM = [float(x) for x in lineSplit[9].split(',')]
                
                simCompIM = [float(x) for x in lineSplit[10].split(',')]
                simContIM = [float(x) for x in lineSplit[11].split(',')]
                
                simCompMS = [float(x) for x in lineSplit[12].split(',')]
                simContMS = [float(x) for x in lineSplit[13].split(',')]
                
                simCompRMS = [float(x) for x in lineSplit[14].split(',')]
                simContRMS = [float(x) for x in lineSplit[15].split(',')]
                
                results[simId] = [bestCompIM, bestContIM, domCompIM, domContIM, simCompIM, simContIM, simCompMS, simContMS, simCompRMS, simContRMS]
                
        print '    Number of test genomes: ' + str(len(genomeIds))
        
        return results
    
    def conditionsPlot(self, results):
        # summarize results for each experimental condition  
        print '  Tabulating results for each experimental condition.'
        itemsProcessed = 0      
        compDataDict = defaultdict(lambda : defaultdict(list))
        contDataDict = defaultdict(lambda : defaultdict(list))
        comps = set()
        conts = set()
        seqLens = set()
        
        compOutliers = defaultdict(list)
        contOutliers = defaultdict(list)
        for simId in results:
            itemsProcessed += 1
            #statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            #sys.stdout.write('%s\r' % statusStr)
            #sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            expCondStr = str(float(comp)) + '-' + str(float(cont)) + '-' + str(int(seqLen))
            
            if genomeId == '647000339' and seqLen == '5000' and comp != '1.00':
                print simId
                print genomeId
                print [('%.1f' % r) for r in results[simId][4]]
                print '*************************'
            
            comps.add(float(comp))
            conts.add(float(cont))
            seqLens.add(int(seqLen))
            
            compDataDict[expCondStr]['best'] += results[simId][0]
            compDataDict[expCondStr]['domain'] += results[simId][2]
            compDataDict[expCondStr]['selected'] += results[simId][4]
            
            for dComp in results[simId][4]:
                compOutliers[expCondStr] += [[dComp, genomeId]]
            
            contDataDict[expCondStr]['best'] += results[simId][1]
            contDataDict[expCondStr]['domain'] += results[simId][3]
            contDataDict[expCondStr]['selected'] += results[simId][5]
            
            for dCont in results[simId][5]:
                contOutliers[expCondStr] += [[dCont, genomeId]]
              
        sys.stdout.write('\n')
        
        print '    There are %d experimental conditions.' % (len(compDataDict))
                
        # plot data
        print '  Plotting results.'
        compData = []
        contData = []
        rowLabels = []
        
        foutComp = open('./experiments/simulation.draft.conditions.comp_outliers.tsv', 'w')
        foutCont = open('./experiments/simulation.draft.conditions.cont_outliers.tsv', 'w')
        for comp in [0.5, 0.7, 0.9]: #sorted(comps):
            for cont in [0.05, 0.1, 0.2]:
                for msStr in ['best', 'selected', 'domain']:
                    for seqLen in [5000]: #sorted(seqLens):
                        rowLabels.append(msStr +': %d%%, %d%%' % (comp*100, cont*100))
                        
                        expCondStr = str(comp) + '-' + str(cont) + '-' + str(seqLen)
                        compData.append(compDataDict[expCondStr][msStr])
                        contData.append(contDataDict[expCondStr][msStr])  
                    
                # report completenes outliers
                foutComp.write(expCondStr)

                compOutliers[expCondStr].sort()
                
                dComps = array([r[0] for r in compOutliers[expCondStr]])
                perc5 = scoreatpercentile(dComps, 5)
                perc95 = scoreatpercentile(dComps, 95)
                print expCondStr, perc5, perc95
                
                outliers = []
                for item in compOutliers[expCondStr]:
                    if item[0] < perc5 or item[0] > perc95:
                        outliers.append(item[1])
                        
                outlierCount = Counter(outliers)
                for genomeId, count in outlierCount.most_common():
                    foutComp.write('\t' + genomeId + ': ' + str(count))
                foutComp.write('\n')
                
                # report contamination outliers
                foutCont.write(expCondStr)

                contOutliers[expCondStr].sort()
                
                dConts = array([r[0] for r in contOutliers[expCondStr]])
                perc5 = scoreatpercentile(dConts, 5)
                perc95 = scoreatpercentile(dConts, 95)
                
                outliers = []
                for item in contOutliers[expCondStr]:
                    if item[0] < perc5 or item[0] > perc95:
                        outliers.append(item[1])
                        
                outlierCount = Counter(outliers)
                for genomeId, count in outlierCount.most_common():
                    foutCont.write('\t' + genomeId + ': ' + str(count))
                foutCont.write('\n')
                
        foutComp.close()
        foutCont.close()
                        
                        
        print 'best:\t%.2f\t%.2f' % (mean(abs(array(compData[0::3]))), mean(abs(array(contData[0::3]))))
        print 'selected:\t%.2f\t%.2f' % (mean(abs(array(compData[1::3]))), mean(abs(array(contData[1::3]))))   
        print 'domain:\t%.2f\t%.2f' % (mean(abs(array(compData[2::3]))), mean(abs(array(contData[2::3]))))   
            
        boxPlot = BoxPlot()
        plotFilename = './experiments/simulation.draft.conditions.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', 'Simulation Conditions', 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 3)
        
    def taxonomicPlots(self, results):
        # summarize results for different taxonomic groups  
        print '  Tabulating results for taxonomic groups.'
        
        metadata = self.img.genomeMetadata()
        
        itemsProcessed = 0      
        compDataDict = defaultdict(lambda : defaultdict(list))
        contDataDict = defaultdict(lambda : defaultdict(list))
        comps = set()
        conts = set()
        seqLens = set()
        
        ranksToProcess = 2
        taxaByRank = [set() for _ in xrange(0, ranksToProcess)]
        
        genomeInTaxon = defaultdict(set)
        for simId in results:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            taxonomy = metadata[genomeId]['taxonomy']
            
            comps.add(float(comp))
            conts.add(float(cont))
            seqLens.add(int(seqLen))
            
            for r in xrange(0, ranksToProcess):
                taxon = taxonomy[r]
                
                if taxon == 'unclassified':
                    continue
                
                taxaByRank[r].add(taxon)
                                
                compDataDict[taxon]['best'] += results[simId][0]
                compDataDict[taxon]['domain'] += results[simId][2]
                compDataDict[taxon]['selected'] += results[simId][4]
                
                contDataDict[taxon]['best'] += results[simId][1]
                contDataDict[taxon]['domain'] += results[simId][3]
                contDataDict[taxon]['selected'] += results[simId][5]
                
                genomeInTaxon[taxon].add(genomeId)
            
        sys.stdout.write('\n')
        
        print '    There are %d taxon.' % (len(compDataDict))
        
        # get list of ordered taxa by rank
        orderedTaxa = []
        for taxa in taxaByRank:
            orderedTaxa += sorted(taxa)
                
        # plot data
        print '  Plotting results.'
        compData = []
        contData = []
        rowLabels = []
        for taxon in orderedTaxa:
            for msStr in ['best', 'selected', 'domain']:
                numGenomes = len(genomeInTaxon[taxon])
                if numGenomes < 10: # skip groups with only a few genomes
                    continue
                
            
                rowLabels.append(msStr + ': ' + taxon + ' (' + str(numGenomes) + ')')
                compData.append(compDataDict[taxon][msStr])
                contData.append(contDataDict[taxon][msStr])        
                
        for i, rowLabel in enumerate(rowLabels):
            print rowLabel + '\t%.2f\t%.2f' % (mean(abs(array(compData[i]))), mean(abs(array(contData[i]))))            
            
        boxPlot = BoxPlot()
        plotFilename = './experiments/simulation.draft.taxonomy.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', None, 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 3)
    
    
    def refinementPlots(self, results):
        # summarize results for different CheckM refinements 
        print '  Tabulating results for different refinements.'
        
        metadata = self.img.genomeMetadata()
        
        itemsProcessed = 0      
        compDataDict = defaultdict(lambda : defaultdict(list))
        contDataDict = defaultdict(lambda : defaultdict(list))
        comps = set()
        conts = set()
        seqLens = set()
        
        ranksToProcess = 2
        taxaByRank = [set() for _ in xrange(0, ranksToProcess)]
        
        genomeInTaxon = defaultdict(set)
        for simId in results:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            taxonomy = metadata[genomeId]['taxonomy']
            
            comps.add(float(comp))
            conts.add(float(cont))
            seqLens.add(int(seqLen))
            
            for r in xrange(0, ranksToProcess):
                taxon = taxonomy[r]
                
                if taxon == 'unclassified':
                    continue
                
                taxaByRank[r].add(taxon)
                
                compDataDict[taxon]['IM'] += results[simId][4]
                compDataDict[taxon]['MS'] += results[simId][6]
                compDataDict[taxon]['RMS'] += results[simId][8]
                
                contDataDict[taxon]['IM'] += results[simId][5]
                contDataDict[taxon]['MS'] += results[simId][7]
                contDataDict[taxon]['RMS'] += results[simId][9]
                
                genomeInTaxon[taxon].add(genomeId)
            
        sys.stdout.write('\n')
        
        print '    There are %d taxon.' % (len(compDataDict))
        
        # get list of ordered taxa by rank
        orderedTaxa = []
        for taxa in taxaByRank:
            orderedTaxa += sorted(taxa)
                
        # plot data
        print '  Plotting results.'
        compData = []
        contData = []
        rowLabels = []
        for taxon in orderedTaxa:
            for refineStr in ['RMS', 'MS', 'IM']:
                numGenomes = len(genomeInTaxon[taxon])
                if numGenomes < 10: # skip groups with only a few genomes
                    continue

                rowLabels.append(refineStr + ': ' + taxon + ' (' + str(numGenomes) + ')')
                compData.append(compDataDict[taxon][refineStr])
                contData.append(contDataDict[taxon][refineStr])       
                
        for i, rowLabel in enumerate(rowLabels):
            print rowLabel + '\t%.2f\t%.2f' % (mean(abs(array(compData[i]))), mean(abs(array(contData[i]))))
            
        boxPlot = BoxPlot()
        plotFilename = './experiments/simulation.draft.refinements.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', None, 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 3)
        
    def run(self):
        # read simulation results
        print '  Reading simulation results.'
        results = self.__readResults('./experiments/simCompare.draft.full.tsv')
                   
        print '\n'         
        self.conditionsPlot(results)
        
        #print '\n'
        #self.taxonomicPlots(results)
        
        #print '\n'
        #self.refinementPlots(results)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simComparePlots = SimComparePlots()
    simComparePlots.run()
