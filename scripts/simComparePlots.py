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
Produce box-and-whisker plots comparing performance of different marker sets (best vs. selected vs. domain).
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
from checkm.util.img import IMG
from checkm.util.taxonomyUtils import rankPrefixes

from numpy import mean, abs, array, std

from scipy.stats import scoreatpercentile

class SimComparePlots(object):
    def __init__(self):
        
        self.plotPrefix = './simulations/simulation.draft.w_refinement_50'
        self.simCompareFile = './simulations/simCompare.draft.w_refinement_50.full.tsv'
        self.simCompareMarkerSetOut = './simulations/simCompare.draft.marker_set_table.w_refinement_50.tsv'
        self.simCompareConditionOut = './simulations/simCompare.draft.condition_table.w_refinement_50.tsv'
        self.simCompareTaxonomyTableOut = './simulations/simCompare.draft.taxonomy_table.w_refinement_50.tsv'
        self.simCompareRefinementTableOut = './simulations/simCompare.draft.refinment_table.w_refinement_50.tsv'
               
        #self.plotPrefix = './simulations/simulation.scaffolds.draft.w_refinement_50'
        #self.simCompareFile = './simulations/simCompare.scaffolds.draft.w_refinement_50.full.tsv'
        #self.simCompareMarkerSetOut = './simulations/simCompare.scaffolds.draft.marker_set_table.w_refinement_50.tsv'
        #self.simCompareConditionOut = './simulations/simCompare.scaffolds.draft.condition_table.w_refinement_50.tsv'
        #self.simCompareTaxonomyTableOut = './simulations/simCompare.scaffolds.draft.taxonomy_table.w_refinement_50.tsv'
        #self.simCompareRefinementTableOut = './simulations/simCompare.scaffolds.draft.refinment_table.w_refinement_50.tsv'
        
        #self.plotPrefix = './simulations/simulation.random_scaffolds.w_refinement_50'
        #self.simCompareFile = './simulations/simCompare.random_scaffolds.w_refinement_50.full.tsv'
        #self.simCompareMarkerSetOut = './simulations/simCompare.random_scaffolds.marker_set_table.w_refinement_50.tsv'
        #self.simCompareConditionOut = './simulations/simCompare.random_scaffolds.condition_table.w_refinement_50.tsv'
        #self.simCompareTaxonomyTableOut = './simulations/simCompare.random_scaffolds.taxonomy_table.w_refinement_50.tsv'
        #self.simCompareRefinementTableOut = './simulations/simCompare.random_scaffolds.refinment_table.w_refinement_50.tsv'
        
        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        
        self.compsToConsider = [0.5, 0.7, 0.8, 0.9] #[0.5, 0.7, 0.8, 0.9]
        self.contsToConsider = [0.05, 0.1, 0.15] #[0.05, 0.1, 0.15]
        
        self.dpi = 1200
  
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
                
                bestCompMS = [float(x) for x in lineSplit[8].split(',')]
                bestContMS = [float(x) for x in lineSplit[9].split(',')]
                                
                domCompIM = [float(x) for x in lineSplit[10].split(',')]
                domContIM = [float(x) for x in lineSplit[11].split(',')]
                
                domCompMS = [float(x) for x in lineSplit[12].split(',')]
                domContMS = [float(x) for x in lineSplit[13].split(',')]
                
                simCompIM = [float(x) for x in lineSplit[14].split(',')]
                simContIM = [float(x) for x in lineSplit[15].split(',')]
                
                simCompMS = [float(x) for x in lineSplit[16].split(',')]
                simContMS = [float(x) for x in lineSplit[17].split(',')]
                
                simCompRMS = [float(x) for x in lineSplit[18].split(',')]
                simContRMS = [float(x) for x in lineSplit[19].split(',')]
                
                results[simId] = [bestCompIM, bestContIM, bestCompMS, bestContMS, domCompIM, domContIM, domCompMS, domContMS, simCompIM, simContIM, simCompMS, simContMS, simCompRMS, simContRMS]
                
        print '    Number of test genomes: ' + str(len(genomeIds))
        
        return results
    
    def markerSets(self, results):
        # summarize results from IM vs MS
        print '  Tabulating results for domain-level marker genes vs marker sets.'
        
        itemsProcessed = 0      
        compDataDict = defaultdict(lambda : defaultdict(list))
        contDataDict = defaultdict(lambda : defaultdict(list))

        genomeIds = set()
        for simId in results:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            genomeIds.add(genomeId)
            expCondStr = str(float(comp)) + '-' + str(float(cont)) + '-' + str(int(seqLen))
            
            compDataDict[expCondStr]['IM'] += results[simId][4]
            compDataDict[expCondStr]['MS'] += results[simId][6]

            contDataDict[expCondStr]['IM'] += results[simId][5]
            contDataDict[expCondStr]['MS'] += results[simId][7]
                
        print '  There are %d unique genomes.' % len(genomeIds)
              
        sys.stdout.write('\n')
        
        print '    There are %d experimental conditions.' % (len(compDataDict))
                
        # plot data
        print '  Plotting results.'
        compData = []
        contData = []
        rowLabels = []
        
        for comp in self.compsToConsider:
            for cont in self.contsToConsider:
                for seqLen in [20000]: 
                    for msStr in ['MS', 'IM']:
                        rowLabels.append(msStr +': %d%%, %d%%' % (comp*100, cont*100))
                        
                        expCondStr = str(comp) + '-' + str(cont) + '-' + str(seqLen)
                        compData.append(compDataDict[expCondStr][msStr])
                        contData.append(contDataDict[expCondStr][msStr])  
                                       
        print 'MS:\t%.2f\t%.2f' % (mean(abs(array(compData[0::2]))), mean(abs(array(contData[0::2]))))
        print 'IM:\t%.2f\t%.2f' % (mean(abs(array(compData[1::2]))), mean(abs(array(contData[1::2]))))   
            
        boxPlot = BoxPlot()
        plotFilename = self.plotPrefix + '.markerSets.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', 'Simulation Conditions', 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 2, dpi = self.dpi)
        
        # print table of results 
        tableOut = open(self.simCompareMarkerSetOut, 'w')
        tableOut.write('Comp. (%)\tCont. (%)\tIM (5kb)\t\tMS (5kb)\t\tIM (20kb)\t\tMS (20kb)\t\tIM (50kb)\t\tMS (50kb)\n')
        
        avgComp = defaultdict(lambda : defaultdict(list))
        avgCont = defaultdict(lambda : defaultdict(list))
        for comp in [0.5, 0.7, 0.8, 0.9, 0.95, 1.0]:
            for cont in [0.0, 0.05, 0.1, 0.15, 0.2]:
                
                tableOut.write('%d\t%d' % (comp*100, cont*100))
                
                for seqLen in [5000, 20000, 50000]:
                    expCondStr = str(comp) + '-' + str(cont) + '-' + str(seqLen)
                     
                    meanCompIM = mean(abs(array(compDataDict[expCondStr]['IM'])))
                    stdCompIM = std(abs(array(compDataDict[expCondStr]['IM'])))
                    meanContIM = mean(abs(array(contDataDict[expCondStr]['IM'])))
                    stdContIM = std(abs(array(contDataDict[expCondStr]['IM'])))
                    
                    avgComp[seqLen]['IM'] += compDataDict[expCondStr]['IM']
                    avgCont[seqLen]['IM'] += contDataDict[expCondStr]['IM']
                    
                    meanCompMS = mean(abs(array(compDataDict[expCondStr]['MS'])))
                    stdCompMS = std(abs(array(compDataDict[expCondStr]['MS'])))
                    meanContMS = mean(abs(array(contDataDict[expCondStr]['MS'])))
                    stdContMS = std(abs(array(contDataDict[expCondStr]['MS'])))
                    
                    avgComp[seqLen]['MS'] += compDataDict[expCondStr]['MS']
                    avgCont[seqLen]['MS'] += contDataDict[expCondStr]['MS']
                    
                    tableOut.write('\t%.1f+/-%.2f\t%.1f+/-%.2f\t%.1f+/-%.2f\t%.1f+/-%.2f' % (meanCompIM, stdCompIM, meanCompMS, stdCompMS, meanContIM, stdContIM, meanContMS, stdContMS))
                tableOut.write('\n')
                
        tableOut.write('\tAverage:')
        for seqLen in [5000, 20000, 50000]: 
            meanCompIM = mean(abs(array(avgComp[seqLen]['IM'])))
            stdCompIM = std(abs(array(avgComp[seqLen]['IM'])))
            meanContIM = mean(abs(array(avgCont[seqLen]['IM'])))
            stdContIM = std(abs(array(avgCont[seqLen]['IM'])))
            
            meanCompMS = mean(abs(array(avgComp[seqLen]['MS'])))
            stdCompMS = std(abs(array(avgComp[seqLen]['MS'])))
            meanContMS = mean(abs(array(avgCont[seqLen]['MS'])))
            stdContMS = std(abs(array(avgCont[seqLen]['MS'])))
            
            tableOut.write('\t%.1f+/-%.2f\t%.1f+/-%.2f\t%.1f+/-%.2f\t%.1f+/-%.2f' % (meanCompIM, stdCompIM, meanCompMS, stdCompMS, meanContIM, stdContIM, meanContMS, stdContMS))
                        
        tableOut.write('\n')     
                
        tableOut.close()
    
    def conditionsPlot(self, results):
        # summarize results for each experimental condition  
        print '  Tabulating results for each experimental condition using marker sets.'
        
        itemsProcessed = 0      
        compDataDict = defaultdict(lambda : defaultdict(list))
        contDataDict = defaultdict(lambda : defaultdict(list))
        comps = set()
        conts = set()
        seqLens = set()
        
        compOutliers = defaultdict(list)
        contOutliers = defaultdict(list)
        
        genomeIds = set()
        for simId in results:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            genomeIds.add(genomeId)
            expCondStr = str(float(comp)) + '-' + str(float(cont)) + '-' + str(int(seqLen))
            
            comps.add(float(comp))
            conts.add(float(cont))
            seqLens.add(int(seqLen))
            
            compDataDict[expCondStr]['best'] += results[simId][2]
            compDataDict[expCondStr]['domain'] += results[simId][6]
            compDataDict[expCondStr]['selected'] += results[simId][10]
            
            for dComp in results[simId][2]:
                compOutliers[expCondStr] += [[dComp, genomeId]]
            
            contDataDict[expCondStr]['best'] += results[simId][3]
            contDataDict[expCondStr]['domain'] += results[simId][7]
            contDataDict[expCondStr]['selected'] += results[simId][11]
            
            for dCont in results[simId][3]:
                contOutliers[expCondStr] += [[dCont, genomeId]]
                
        print '  There are %d unique genomes.' % len(genomeIds)
              
        sys.stdout.write('\n')
        
        print '    There are %d experimental conditions.' % (len(compDataDict))
                
        # plot data
        print '  Plotting results.'
        compData = []
        contData = []
        rowLabels = []
        
        foutComp = open('./simulations/simulation.scaffolds.draft.comp_outliers.domain.tsv', 'w')
        foutCont = open('./simulations/simulation.scaffolds.draft.cont_outliers.domain.tsv', 'w')
        for comp in self.compsToConsider:
            for cont in self.contsToConsider:
                for msStr in ['best', 'selected', 'domain']:
                    for seqLen in [20000]: 
                        rowLabels.append(msStr +': %d%%, %d%%' % (comp*100, cont*100))
                        
                        expCondStr = str(comp) + '-' + str(cont) + '-' + str(seqLen)
                        compData.append(compDataDict[expCondStr][msStr])
                        contData.append(contDataDict[expCondStr][msStr])  
                    
                # report completenes outliers
                foutComp.write(expCondStr)

                compOutliers[expCondStr].sort()
                
                dComps = array([r[0] for r in compOutliers[expCondStr]])
                perc1 = scoreatpercentile(dComps, 1)
                perc99 = scoreatpercentile(dComps, 99)
                print expCondStr, perc1, perc99
                
                foutComp.write('\t%.2f\t%.2f' % (perc1, perc99))
                
                outliers = []
                for item in compOutliers[expCondStr]:
                    if item[0] < perc1 or item[0] > perc99:
                        outliers.append(item[1])
                        
                outlierCount = Counter(outliers)
                for genomeId, count in outlierCount.most_common():
                    foutComp.write('\t' + genomeId + ': ' + str(count))
                foutComp.write('\n')
                
                # report contamination outliers
                foutCont.write(expCondStr)

                contOutliers[expCondStr].sort()
                
                dConts = array([r[0] for r in contOutliers[expCondStr]])
                perc1 = scoreatpercentile(dConts, 1)
                perc99 = scoreatpercentile(dConts, 99)
                
                foutCont.write('\t%.2f\t%.2f' % (perc1, perc99))
                
                outliers = []
                for item in contOutliers[expCondStr]:
                    if item[0] < perc1 or item[0] > perc99:
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
        plotFilename = self.plotPrefix + '.conditions.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', 'Simulation Conditions', 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 3, dpi = self.dpi)
        
        
        # print table of results 
        tableOut = open(self.simCompareConditionOut, 'w')
        tableOut.write('Comp. (%)\tCont. (%)\tbest (5kb)\t\tselected (5kb)\t\tdomain (5kb)\t\tbest (20kb)\t\tselected (20kb)\t\tdomain (20kb)\t\tbest (50kb)\t\tselected (50kb)\t\tdomain (50kb)\n')
        
        avgComp = defaultdict(lambda : defaultdict(list))
        avgCont = defaultdict(lambda : defaultdict(list))
        for comp in [0.5, 0.7, 0.8, 0.9, 0.95, 1.0]:
            for cont in [0.0, 0.05, 0.1, 0.15, 0.2]:
                
                tableOut.write('%d\t%d' % (comp*100, cont*100))
                
                for seqLen in [5000, 20000, 50000]:
                    expCondStr = str(comp) + '-' + str(cont) + '-' + str(seqLen)
                   
                    meanCompD = mean(abs(array(compDataDict[expCondStr]['domain'])))
                    stdCompD = std(abs(array(compDataDict[expCondStr]['domain'])))
                    meanContD = mean(abs(array(contDataDict[expCondStr]['domain'])))
                    stdContD = std(abs(array(contDataDict[expCondStr]['domain'])))
                    
                    avgComp[seqLen]['domain'] += compDataDict[expCondStr]['domain']
                    avgCont[seqLen]['domain'] += contDataDict[expCondStr]['domain']
                    
                    meanCompS = mean(abs(array(compDataDict[expCondStr]['selected'])))
                    stdCompS = std(abs(array(compDataDict[expCondStr]['selected'])))
                    meanContS = mean(abs(array(contDataDict[expCondStr]['selected'])))
                    stdContS = std(abs(array(contDataDict[expCondStr]['selected'])))
                    
                    avgComp[seqLen]['selected'] += compDataDict[expCondStr]['selected']
                    avgCont[seqLen]['selected'] += contDataDict[expCondStr]['selected']
                    
                    meanCompB = mean(abs(array(compDataDict[expCondStr]['best'])))
                    stdCompB = std(abs(array(compDataDict[expCondStr]['best'])))
                    meanContB = mean(abs(array(contDataDict[expCondStr]['best'])))
                    stdContB = std(abs(array(contDataDict[expCondStr]['best'])))
                    
                    avgComp[seqLen]['best'] += compDataDict[expCondStr]['best']
                    avgCont[seqLen]['best'] += contDataDict[expCondStr]['best']
                    
                    tableOut.write('\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f' % (meanCompD, meanCompS, meanCompB, meanContD, meanContS, meanContB))
                tableOut.write('\n')
                
        tableOut.write('\tAverage:')
        for seqLen in [5000, 20000, 50000]: 
            meanCompD = mean(abs(array(avgComp[seqLen]['domain'])))
            stdCompD = std(abs(array(avgComp[seqLen]['domain'])))
            meanContD = mean(abs(array(avgCont[seqLen]['domain'])))
            stdContD = std(abs(array(avgCont[seqLen]['domain'])))
            
            meanCompS = mean(abs(array(avgComp[seqLen]['selected'])))
            stdCompS = std(abs(array(avgComp[seqLen]['selected'])))
            meanContS = mean(abs(array(avgCont[seqLen]['selected'])))
            stdContS = std(abs(array(avgCont[seqLen]['selected'])))
            
            meanCompB = mean(abs(array(avgComp[seqLen]['best'])))
            stdCompB = std(abs(array(avgComp[seqLen]['best'])))
            meanContB = mean(abs(array(avgCont[seqLen]['best'])))
            stdContB = std(abs(array(avgCont[seqLen]['best'])))
            
            tableOut.write('\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f' % (meanCompD, meanCompS, meanCompB, meanContD, meanContS, meanContB))
                        
        tableOut.write('\n')     
                
        tableOut.close()
        
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
        
        ranksToProcess = 3
        taxaByRank = [set() for _ in xrange(0, ranksToProcess)]
        
        overallComp = []
        overallCont = []
                
        genomeInTaxon = defaultdict(set)
        testCases = 0
        for simId in results:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            
            if seqLen != '20000':
                continue
            
            if str(float(comp)) in ['0.5', '0.7', '0.8', '0.9'] and str(float(cont)) in ['0.05', '0.10', '0.1', '0.15']:
                print comp, cont
                taxonomy = metadata[genomeId]['taxonomy']
                
                testCases += 1
                
                comps.add(float(comp))
                conts.add(float(cont))
                seqLens.add(int(seqLen))
                
                overallComp += results[simId][10]
                overallCont += results[simId][11]
                
                for r in xrange(0, ranksToProcess):
                    taxon = taxonomy[r]
                    
                    if r == 0 and taxon == 'unclassified':
                        print '*****************************Unclassified at domain-level*****************'
                        continue
                    
                    if taxon == 'unclassified':
                        continue
                    
                    taxon = rankPrefixes[r] + taxon
                    
                    taxaByRank[r].add(taxon)
                                                    
                    compDataDict[taxon]['best'] += results[simId][2]
                    compDataDict[taxon]['domain'] += results[simId][6]
                    compDataDict[taxon]['selected'] += results[simId][10]
                    
                    contDataDict[taxon]['best'] += results[simId][3]
                    contDataDict[taxon]['domain'] += results[simId][7]
                    contDataDict[taxon]['selected'] += results[simId][11]
                    
                    genomeInTaxon[taxon].add(genomeId)
            
        sys.stdout.write('\n')
        
        print 'Test cases', testCases
        
        print ''        
        print 'Creating plots for:'
        print '  comps = ', comps
        print '  conts = ', conts
        
        print ''
        print '    There are %d taxa.' % (len(compDataDict))
        
        print ''
        print '  Overall bias:'
        print '    Selected comp: %.2f' % mean(overallComp)
        print '    Selected cont: %.2f' % mean(overallCont)
        
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
                  
        # print taxonomic table of results organized by class
        taxonomyTableOut = open(self.simCompareTaxonomyTableOut, 'w')
        for taxon in orderedTaxa:
            numGenomes = len(genomeInTaxon[taxon])
            if numGenomes < 2: # skip groups with only a few genomes
                continue
                
            taxonomyTableOut.write(taxon + '\t' + str(numGenomes))
            for msStr in ['domain', 'selected']:                
                meanTaxonComp = mean(abs(array(compDataDict[taxon][msStr])))
                stdTaxonComp = std(abs(array(compDataDict[taxon][msStr])))
                meanTaxonCont = mean(abs(array(contDataDict[taxon][msStr])))
                stdTaxonCont = std(abs(array(contDataDict[taxon][msStr])))
                
                taxonomyTableOut.write('\t%.1f +/- %.2f\t%.1f +/- %.2f' % (meanTaxonComp, stdTaxonComp, meanTaxonCont, stdTaxonCont))
            taxonomyTableOut.write('\n')
        taxonomyTableOut.close()
        
        # create box plot
        boxPlot = BoxPlot()
        plotFilename = self.plotPrefix +  '.taxonomy.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', None, 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 3, dpi = self.dpi)
    
    
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
        
        ranksToProcess = 3
        taxaByRank = [set() for _ in xrange(0, ranksToProcess)]
        
        overallCompIM = []
        overallContIM = [] 
        
        overallCompMS = []
        overallContMS = [] 
        
        overallCompRMS = []
        overallContRMS = [] 
        
        genomeInTaxon = defaultdict(set)
        
        testCases = 0
        for simId in results:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, seqLen, comp, cont = simId.split('-')
            taxonomy = metadata[genomeId]['taxonomy']
            
            if float(comp) < 0.7 or float(cont) > 0.1:
                continue
            
            comps.add(float(comp))
            conts.add(float(cont))
            seqLens.add(int(seqLen))
            
            overallCompIM.append(results[simId][8])
            overallContIM.append(results[simId][9])
            
            overallCompMS.append(results[simId][10])
            overallContMS.append(results[simId][11])
            
            overallCompRMS.append(results[simId][12])
            overallContRMS.append(results[simId][13])
            
            for r in xrange(0, ranksToProcess):
                taxon = taxonomy[r]
                
                if taxon == 'unclassified':
                    continue
                
                taxaByRank[r].add(taxon)
                
                compDataDict[taxon]['IM'] += results[simId][8]
                compDataDict[taxon]['MS'] += results[simId][10]
                compDataDict[taxon]['RMS'] += results[simId][12]
                
                contDataDict[taxon]['IM'] += results[simId][9]
                contDataDict[taxon]['MS'] += results[simId][11]
                contDataDict[taxon]['RMS'] += results[simId][13]
                                
                genomeInTaxon[taxon].add(genomeId)
            
        sys.stdout.write('\n')
        
        print 'Creating plots for:'
        print '  comps = ', comps
        print '  conts = ', conts
        
        print ''
        print '    There are %d taxon.' % (len(compDataDict))
        print ''
        print 'Percentage change MS-IM comp: %.4f' % ((mean(abs(array(overallCompMS))) - mean(abs(array(overallCompIM)))) * 100 / mean(abs(array(overallCompIM))))
        print 'Percentage change MS-IM cont: %.4f' % ((mean(abs(array(overallContMS))) - mean(abs(array(overallContIM)))) * 100 / mean(abs(array(overallContIM))))
        print ''
        print 'Percentage change RMS-MS comp: %.4f' % ((mean(abs(array(overallCompRMS))) - mean(abs(array(overallCompMS)))) * 100 / mean(abs(array(overallCompIM))))
        print 'Percentage change RMS-MS cont: %.4f' % ((mean(abs(array(overallContRMS))) - mean(abs(array(overallContMS)))) * 100 / mean(abs(array(overallContIM))))
        
        print ''
        
        # get list of ordered taxa by rank
        orderedTaxa = []
        for taxa in taxaByRank:
            orderedTaxa += sorted(taxa)
             
        # print table of results organized by class
        refinmentTableOut = open(self.simCompareRefinementTableOut, 'w')
        for taxon in orderedTaxa:
            numGenomes = len(genomeInTaxon[taxon])
            if numGenomes < 2: # skip groups with only a few genomes
                continue
                
            refinmentTableOut.write(taxon + '\t' + str(numGenomes))
            for refineStr in ['IM', 'MS']:               
                meanTaxonComp = mean(abs(array(compDataDict[taxon][refineStr])))
                stdTaxonComp = std(abs(array(compDataDict[taxon][refineStr])))
                meanTaxonCont = mean(abs(array(contDataDict[taxon][refineStr])))
                stdTaxonCont = std(abs(array(contDataDict[taxon][refineStr])))
                
                refinmentTableOut.write('\t%.1f +/- %.2f\t%.1f +/- %.2f' % (meanTaxonComp, stdTaxonComp, meanTaxonCont, stdTaxonCont))
            
            perCompChange = (mean(abs(array(compDataDict[taxon]['IM']))) - meanTaxonComp) * 100 / mean(abs(array(compDataDict[taxon]['IM'])))
            perContChange = (mean(abs(array(contDataDict[taxon]['IM']))) - meanTaxonCont) * 100 / mean(abs(array(contDataDict[taxon]['IM'])))
            refinmentTableOut.write('\t%.2f\t%.2f\n' % (perCompChange, perContChange))
        refinmentTableOut.close()
       
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
        plotFilename = self.plotPrefix + '.refinements.png'
        boxPlot.plot(plotFilename, compData, contData, rowLabels, 
                        r'$\Delta$' + ' % Completion', None, 
                        r'$\Delta$' + ' % Contamination', None,
                        rowsPerCategory = 3, dpi = self.dpi)
        
    def run(self):
        # read simulation results
        print '  Reading simulation results.'
        results = self.__readResults(self.simCompareFile)
        
        print '\n'         
        #self.markerSets(results)
                   
        print '\n'         
        #self.conditionsPlot(results)
        
        #print '\n'
        self.taxonomicPlots(results)
        
        print '\n'
        #self.refinementPlots(results)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simComparePlots = SimComparePlots()
    simComparePlots.run()
