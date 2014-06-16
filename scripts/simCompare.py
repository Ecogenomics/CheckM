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
Compare best performing marker set to marker set selected via simulation.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import gzip

from collections import defaultdict

from numpy import array, mean, std, percentile

import dendropy
from  dendropy.dataobject.taxon import Taxon

class SimCompare(object):
    def __init__(self):
        self.resultsSummaryFile = './simulations/simulation.draft.summary.tsv'
        self.resultsFullFile = './simulations/simulation.draft.tsv.gz'
        self.simCompareSummaryOut = './simulations/simCompare.draft.summary.tsv'
        self.simCompareFullOut = './simulations/simCompare.draft.full.tsv'
        
        #self.resultsSummaryFile = './simulations/simulation.scaffolds.draft.summary.tsv'
        #self.resultsFullFile = './simulations/simulation.scaffolds.draft.tsv.gz'
        #self.simCompareSummaryOut = './simulations/simCompare.scaffolds.draft.summary.tsv'
        #self.simCompareFullOut = './simulations/simCompare.scaffolds.draft.full.tsv'
            
    def __bestMarkerSet(self, simId, simResults):
        """Get stats for best marker set."""
        bestUID, numDescendantsBest, dCompBest, dContBest = self.__domainMarkerSet(simId, simResults)
        curBest = dCompBest + dContBest
        for uid, results in simResults[simId].iteritems():
            numDescendants, dComp, dCont, _, _, _, _, _, _ = results
            if numDescendants < 10:
                continue 
            
            if (dComp + dCont) < curBest:
                numDescendantsBest = numDescendants
                dCompBest = dComp
                dContBest = dCont
                bestUID = uid
                curBest = dComp + dCont
                
        return bestUID, numDescendantsBest, dCompBest, dContBest
    
    def __domainMarkerSet(self, simId, simResults):
        if 'UID2' in simResults[simId]:
            return 'UID2', simResults[simId]['UID2'][0], simResults[simId]['UID2'][1],  simResults[simId]['UID2'][2] # archaea
        
        return 'UID203', simResults[simId]['UID203'][0], simResults[simId]['UID203'][1],  simResults[simId]['UID203'][2] # bacteria
    
    def __inferredMarkerSet(self, tree, genomeId, inferredMarkerSet):
        curNode = tree.find_node_with_taxon_label('IMG_' + genomeId)
        
        curNode = curNode.parent_node
        while curNode != None:
            uniqueId = curNode.label.split('|')[0] 
            if inferredMarkerSet.get(uniqueId, 'NA') != 'NA':
                return inferredMarkerSet[uniqueId]
            
            # return domain-specific set if reached that far
            if uniqueId == 'UID2':
                return 'UID2'
            elif uniqueId == 'UID203':
                return 'UID203'
            
            curNode = curNode.parent_node
            
        return None
    
    def __readSummaryResults(self, filename):
        summaryResults = defaultdict(dict)
        genomeIds = set()
        with open(filename) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                simId = lineSplit[0] + '-' + lineSplit[1] + '-' + lineSplit[2] + '-' + lineSplit[3]
                genomeIds.add(lineSplit[0])
                uid = lineSplit[5].split('|')[0].strip()
                numDescendants = int(lineSplit[6])
                
                compUnmodified = float(lineSplit[7].rstrip())
                contUnmodified = float(lineSplit[8].rstrip())
                
                compIM = float(lineSplit[9].rstrip())
                contIM = float(lineSplit[11].rstrip())
                
                compMS = float(lineSplit[13].rstrip())
                contMS = float(lineSplit[15].rstrip())
                
                compRMS = float(lineSplit[21].rstrip())
                contRMS = float(lineSplit[23].rstrip())
                
                summaryResults[simId][uid] = [numDescendants, compIM, contIM, compMS, contMS, compRMS, contRMS, compUnmodified, contUnmodified]
                
        print '    Number of test genomes: ' + str(len(genomeIds))
        
        return summaryResults
    
    def __readFullResults(self, filename):
        summaryResults = defaultdict(dict)
        genomeIds = set()
        with gzip.open(filename) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                simId = lineSplit[0] + '-' + lineSplit[1] + '-' + lineSplit[2] + '-' + lineSplit[3]
                genomeIds.add(lineSplit[0])
                uid = lineSplit[5].split('|')[0].strip()
                numDescendants = int(lineSplit[6])
                
                compIM = lineSplit[9].rstrip()
                contIM = lineSplit[10].rstrip()
                
                compMS = lineSplit[11].rstrip()
                contMS = lineSplit[12].rstrip()
                
                compRMS = lineSplit[15].rstrip()
                contRMS = lineSplit[16].rstrip()
            
                summaryResults[simId][uid] = [numDescendants, compIM, contIM, compMS, contMS, compRMS, contRMS]
                
        print '    Number of test genomes: ' + str(len(genomeIds))
        
        return summaryResults

    def run(self):
        print '\n  Reading reference genome tree.'
        treeFile = os.path.join('/srv/whitlam/bio/db/checkm/genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # read simulation results        
        print '  Reading summary simulation results.'
        summaryResults = self.__readSummaryResults(self.resultsSummaryFile)
        print '  Reading full simulation results.'
        fullResults = self.__readFullResults(self.resultsFullFile)
                
        # reading marker set inferred via simulation
        inferredMarkerSet = {}
        with open('./simulations/simInferBestMarkerSet.tsv') as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                uid = lineSplit[0].split('|')[0].strip()
                markerSetId = lineSplit[1].split('|')[0].strip()
                inferredMarkerSet[uid] = markerSetId
                
        # determine best marker set
        print '  Comparing best marker set to selected marker set.'
        foutSummary = open(self.simCompareSummaryOut, 'w')
        foutSummary.write('Sim Id\tBest Id\tDescendants\tSelected Id\tDescendants\tDomain Id')
        foutSummary.write('\tdComp Best\tdContBest\tComp domain\tCont domain\tComp sim\tCont sim\tComp sim [MS]\tCont sim [MS]\tComp sim [RMS]\tCont sim [RMS]\n')
        
        foutFull = open(self.simCompareFullOut, 'w')
        foutFull.write('Sim Id\tBest Id\tDescendants\tSelected Id\tDescendants\tDomain Id\tdComp Best [IM]\tdCont Best [IM]\tdComp domain [IM]\tdCont domain [IM]\tdComp sim [IM]\tdCont sim [IM]\tdComp sim [MS]\tdCont sim [MS]\tdComp sim [RMS]\tdCont sim [RMS]\n')
        
        itemsProcessed = 0
        dCompSimBestList = defaultdict(lambda : defaultdict(list))
        dContSimBestList = defaultdict(lambda : defaultdict(list))
        
        dCompDomList = defaultdict(lambda : defaultdict(list))
        dContDomList = defaultdict(lambda : defaultdict(list))
        
        dCompSimList = defaultdict(lambda : defaultdict(list))
        dContSimList = defaultdict(lambda : defaultdict(list))
        
        
        # DEBUG:
        unmodifiedComp = []
        unmodifiedCont = []
        incompleteGenomesSelected = set()
        contaminatedGenomesSelected = set()
        incompleteGenomesDomain = set()
        contaminatedGenomesDomain = set()
        incompleteGenomesBest = set()
        contaminatedGenomesBest = set()
        totalGenomes = set()
        
        
        dCompDomOverall = []
        dContDomOverall = []
        dCompSelectedOverall = []
        dContSelectedOverall = []
                  
        domBetter = 0
        simBetter = 0
        
        msBetter = 0
        rmsBetter = 0
        
        imDomBetter = 0
        msDomBetter = 0
        
        briefSummaryOut = open('./simulations/briefSummaryOut.tsv', 'w')
        for simId in summaryResults:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(summaryResults), float(itemsProcessed)*100/len(summaryResults))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, _, comp, cont = simId.split('-')
            
            #if genomeId != '2529293088':
            #    continue

            # get results for best performing marker set
            bestUID, numDescendantsBest, dCompBestSummary, dContBestSummary = self.__bestMarkerSet(simId, summaryResults)
            _, dCompBestFull, dContBestFull, _, _, _, _ = fullResults[simId][bestUID]
                        
            # get results for domain-level marker set
            domainUID, _, dCompDomSummary, dContDomSummary = self.__domainMarkerSet(simId, summaryResults)
            _, dCompDomFull, dContDomFull, dCompDomFullMS, dContDomFullMS, _, _ = fullResults[simId][domainUID]
            
            # get results for selected marker set
            simUID = self.__inferredMarkerSet(tree, genomeId, inferredMarkerSet)
            numDescendantsSim, dCompSimSummaryIM, dContSimSummaryIM, dCompSimSummaryMS, dContSimSummaryMS, dCompSimSummaryRMS, dContSimSummaryRMS, _, _ = summaryResults[simId][simUID]
            _, dCompSimFullIM, dContSimFullIM, dCompSimFullMS, dContSimFullMS, dCompSimFullRMS, dContSimFullRMS = fullResults[simId][simUID]
               
            for a, b, c, d, e, f in zip(dCompDomFullMS.split(','), dCompSimFullMS.split(','), dCompSimFullRMS.split(','), dContDomFullMS.split(','), dContSimFullMS.split(','), dContSimFullRMS.split(',')):
                briefSummaryOut.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (float(a), float(b), float(c), float(d), float(e), float(f)))

            for a, b in zip(dCompDomFull.split(','), dCompDomFullMS.split(',')):
                if abs(abs(float(a)) - abs(float(b))) > 0.1:             
                    if abs(float(a)) < abs(float(b)):
                        imDomBetter += 1
                    else:
                        msDomBetter += 1
            
            for a, b in zip(dCompDomFullMS.split(','), dCompSimFullMS.split(',')):
                if abs(abs(float(a)) - abs(float(b))) > 0.1:             
                    if abs(float(a)) < abs(float(b)):
                        domBetter += 1
                    else:
                        simBetter += 1
              
            for a, b in zip(dCompSimFullMS.split(','), dCompSimFullRMS.split(',')):
                if abs(abs(float(a)) - abs(float(b))) > 0.1:    
                    if abs(float(a)) < abs(float(b)):
                        rmsBetter += 1
                    else:
                        msBetter += 1
            
            # DEBUG
            unmodifiedComp.append(summaryResults[simId][bestUID][7])
            unmodifiedCont.append(summaryResults[simId][bestUID][8])
            
            totalGenomes.add(genomeId)
            if summaryResults[simId][simUID][7] < 95:
                incompleteGenomesSelected.add(genomeId)
            if summaryResults[simId][simUID][8] > 5:
                contaminatedGenomesSelected.add(genomeId)
                 
            if summaryResults[simId][domainUID][7] < 95:
                incompleteGenomesDomain.add(genomeId)
            if summaryResults[simId][domainUID][8] > 5:
                contaminatedGenomesDomain.add(genomeId)
                
            if summaryResults[simId][bestUID][7] < 95:
                incompleteGenomesBest.add(genomeId)
            if summaryResults[simId][bestUID][8] > 5:
                contaminatedGenomesBest.add(genomeId)    
                
                
                
            dCompDomOverall.append(dCompDomSummary)
            dContDomOverall.append(dContDomSummary)
            dCompSelectedOverall.append(dCompSimSummaryIM)
            dContSelectedOverall.append(dContSimSummaryIM)
                
                   
            
            dCompSimBestList[comp][cont].append(dCompSimSummaryIM - dCompBestSummary)
            dContSimBestList[comp][cont].append(dContSimSummaryIM - dContBestSummary)
            dCompDomList[comp][cont].append(dCompDomSummary)
            dContDomList[comp][cont].append(dContDomSummary)
            dCompSimList[comp][cont].append(dCompSimSummaryIM)
            dContSimList[comp][cont].append(dContSimSummaryIM)
                
            foutSummary.write('%s\t%s\t%d\t%s\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (simId, bestUID, numDescendantsBest, simUID, numDescendantsSim, domainUID, 
                                                                                                        dCompSimSummaryIM - dCompBestSummary, dContSimSummaryIM - dContBestSummary, 
                                                                                                        dCompDomSummary, dContDomSummary, 
                                                                                                        dCompSimSummaryIM, dContSimSummaryIM, dCompSimSummaryMS, dContSimSummaryMS, dCompSimSummaryRMS, dContSimSummaryRMS))
                                
            foutFull.write('%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (simId, bestUID, numDescendantsBest, simUID, numDescendantsSim, domainUID, 
                                                                                                 dCompBestFull, dContBestFull, 
                                                                                                 dCompDomFull, dContDomFull, 
                                                                                                 dCompSimFullIM, dContSimFullIM, dCompSimFullMS, dContSimFullMS, dCompSimFullRMS, dContSimFullRMS))

        foutSummary.close()
        foutFull.close()
                
        sys.stdout.write('\n')
        
        print '\nOverall results:'
        for comp in sorted(dCompSimBestList.keys()):
            for cont in sorted(dCompSimBestList[comp].keys()):
                print 'Comp: %s, cont: %s' % (comp, cont)
                
                dCompSimBest = array(dCompSimBestList[comp][cont])
                dContSimBest = array(dContSimBestList[comp][cont])
                dCompDom = array(dCompDomList[comp][cont])
                dContDom = array(dContDomList[comp][cont])
                dCompSim = array(dCompSimList[comp][cont])
                dContSim = array(dContSimList[comp][cont])
                
                print 'delta Comp Sim - Best: %.3f +/- %.3f, %.3f, %.3f' % (mean(abs(dCompSimBest)), std(abs(dCompSimBest)), percentile(dCompSimBest, 10), percentile(dCompSimBest, 90))
                print 'delta Cont Sim - Best: %.3f +/- %.3f, %.3f, %.3f' % (mean(abs(dContSimBest)), std(abs(dContSimBest)), percentile(dContSimBest, 10), percentile(dContSimBest, 90))
                print 'Comp Dom: %.3f +/- %.3f, %.3f, %.3f' % (mean(abs(dCompDom)), std(abs(dCompDom)), percentile(dCompDom, 10), percentile(dCompDom, 90))
                print 'Cont Dom: %.3f +/- %.3f, %.3f, %.3f' % (mean(abs(dContDom)), std(abs(dContDom)), percentile(dContDom, 10), percentile(dContDom, 90))
                print 'Comp Sim: %.3f +/- %.3f, %.3f, %.3f' % (mean(abs(dCompSim)), std(abs(dCompSim)), percentile(dCompSim, 10), percentile(dCompSim, 90))
                print 'Cont Sim: %.3f +/- %.3f, %.3f, %.3f' % (mean(abs(dContSim)), std(abs(dContSim)), percentile(dContSim, 10), percentile(dContSim, 90))
                print '\n'
                
        print 'Unmodified comp for best ms: mean = %.1f, 5th = %.1f, 95th = %.1f, std = %.2f, min = %.1f' % (mean(unmodifiedComp), percentile(unmodifiedComp, 5), percentile(unmodifiedComp, 95), std(unmodifiedComp), min(unmodifiedComp))
        print 'Unmodified cont for best ms: mean = %.1f, 5th = %.1f, 95th = %.1f, std = %.2f, max = %.1f' % (mean(unmodifiedCont), percentile(unmodifiedCont, 5), percentile(unmodifiedCont, 95), std(unmodifiedCont), max(unmodifiedCont))
        
        print 'There are %d of %d (%.2f%%) genomes with a comp < 95%% on the domain ms.' % (len(incompleteGenomesDomain), len(totalGenomes), len(incompleteGenomesDomain)*100.0/len(totalGenomes))
        print 'There are %d of %d (%.2f%%) genomes with a cont > 5%% on the domain ms.' % (len(contaminatedGenomesDomain), len(totalGenomes), len(contaminatedGenomesDomain)*100.0/len(totalGenomes))
        
        print 'There are %d of %d (%.2f%%) genomes with a comp < 95%% on the selected ms.' % (len(incompleteGenomesSelected), len(totalGenomes), len(incompleteGenomesSelected)*100.0/len(totalGenomes))
        print 'There are %d of %d (%.2f%%) genomes with a cont > 5%% on the selected ms.' % (len(contaminatedGenomesSelected), len(totalGenomes), len(contaminatedGenomesSelected)*100.0/len(totalGenomes))
        
        print 'There are %d of %d (%.2f%%) genomes with a comp < 95%% on the best ms.' % (len(incompleteGenomesBest), len(totalGenomes), len(incompleteGenomesBest)*100.0/len(totalGenomes))
        print 'There are %d of %d (%.2f%%) genomes with a cont > 5%% on the best ms.' % (len(contaminatedGenomesBest), len(totalGenomes), len(contaminatedGenomesBest)*100.0/len(totalGenomes))
                
        #for genomeId in contaminatedGenomesDomain:
        #    os.system('cp /srv/db/img/07042014/genomes/' + genomeId + '/' + genomeId + '.fna ./genome_test/incomplete_dom10')
        
        print 'Completeness for domain and selected: %.2f +/- %.2f, %.2f +/- %.2f' % (mean(abs(array(dCompDomOverall))), std(abs(array(dCompDomOverall))), mean(abs(array(dCompSelectedOverall))), std(abs(array(dCompSelectedOverall))))
        print 'Contamination for domain and selected: %.2f +/- %.2f, %.2f +/- %.2f' % (mean(abs(array(dContDomOverall))), std(abs(array(dContDomOverall))), mean(abs(array(dContSelectedOverall))), std(abs(array(dContSelectedOverall))))
            
            
        print ''
        print domBetter, simBetter
        print 'Domain better: %.2f' % (float(domBetter)*100/(domBetter+simBetter))
        print 'Sim better: %.2f' % (float(simBetter)*100/(domBetter+simBetter))
        
        print ''
        print msBetter, rmsBetter
        print 'MS better: %.2f' % (float(msBetter)*100/(msBetter+rmsBetter))
        print 'RMS better: %.2f' % (float(rmsBetter)*100/(msBetter+rmsBetter))
        
        print ''
        print imDomBetter, msDomBetter
        print 'IM Domain better: %.2f' % (float(imDomBetter)*100/(imDomBetter+msDomBetter))
        print 'MS Domain better: %.2f' % (float(msDomBetter)*100/(imDomBetter+msDomBetter))
        
        briefSummaryOut.close()
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simCompare = SimCompare()
    simCompare.run()
