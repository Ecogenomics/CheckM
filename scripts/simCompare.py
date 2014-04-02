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
        pass
    
    def __bestMarkerSet(self, simId, simResults):
        """Get stats for best marker set."""
        bestUID, numDescendantsBest, dCompBest, dContBest = self.__domainMarkerSet(simId, simResults)
        curBest = dCompBest + dContBest
        for uid, results in simResults[simId].iteritems():
            numDescendants, dComp, dCont, _, _, _, _ = results
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
        
        return 'UID170', simResults[simId]['UID170'][0], simResults[simId]['UID170'][1],  simResults[simId]['UID170'][2] # bacteria
    
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
            elif uniqueId == 'UID170':
                return 'UID170'
            
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
                
                compIM = float(lineSplit[9].rstrip())
                contIM = float(lineSplit[11].rstrip())
                
                compMS = float(lineSplit[13].rstrip())
                contMS = float(lineSplit[15].rstrip())
                
                compRMS = float(lineSplit[21].rstrip())
                contRMS = float(lineSplit[23].rstrip())
                
                summaryResults[simId][uid] = [numDescendants, compIM, contIM, compMS, contMS, compRMS, contRMS]
                
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
        treeFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # read simulation results
        print '  Reading summary simulation results.'
        summaryResults = self.__readSummaryResults('./experiments/simulation.draft.summary.tsv')
        print '  Reading full simulation results.'
        fullResults = self.__readFullResults('./experiments/simulation.draft.tsv.gz')
                
        # reading marker set inferred via simulation
        inferredMarkerSet = {}
        with open('./experiments/simInferBestMarkerSet.tsv') as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                uid = lineSplit[0].split('|')[0].strip()
                markerSetId = lineSplit[1].split('|')[0].strip()
                inferredMarkerSet[uid] = markerSetId
                
        # determine best marker set
        print '  Comparing best marker set to selected marker set.'
        foutSummary = open('./experiments/simCompare.draft.summary.tsv', 'w')
        foutSummary.write('Sim Id\tBest Id\tDescendants\tSelected Id\tDescendants\tDomain Id')
        foutSummary.write('\tdComp Best\tdContBest\tComp domain\tCont domain\tComp sim\tCont sim\tComp sim [MS]\tCont sim [MS]\tComp sim [RMS]\tCont sim [RMS]\n')
        
        foutFull = open('./experiments/simCompare.draft.full.tsv', 'w')
        foutFull.write('Sim Id\tBest Id\tDescendants\tSelected Id\tDescendants\tDomain Id\tdComp Best [IM]\tdCont Best [IM]\tdComp domain [IM]\tdCont domain [IM]\tdComp sim [IM]\tdCont sim [IM]\tdComp sim [MS]\tdCont sim [MS]\tdComp sim [RMS]\tdCont sim [RMS]\n')
        
        itemsProcessed = 0
        dCompSimBestList = defaultdict(lambda : defaultdict(list))
        dContSimBestList = defaultdict(lambda : defaultdict(list))
        
        dCompDomList = defaultdict(lambda : defaultdict(list))
        dContDomList = defaultdict(lambda : defaultdict(list))
        
        dCompSimList = defaultdict(lambda : defaultdict(list))
        dContSimList = defaultdict(lambda : defaultdict(list))
        
        for simId in summaryResults:
            itemsProcessed += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(summaryResults), float(itemsProcessed)*100/len(summaryResults))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, _, comp, cont = simId.split('-')

            # get results for best performing marker set
            bestUID, numDescendantsBest, dCompBestSummary, dContBestSummary = self.__bestMarkerSet(simId, summaryResults)
            _, dCompBestFull, dContBestFull, _, _, _, _ = fullResults[simId][bestUID]
            
            # get results for domain-level marker set
            domainUID, _, dCompDomSummary, dContDomSummary = self.__domainMarkerSet(simId, summaryResults)
            _, dCompDomFull, dContDomFull, _, _, _, _ = fullResults[simId][domainUID]
            
            # get results for selected marker set
            simUID = self.__inferredMarkerSet(tree, genomeId, inferredMarkerSet)
            numDescendantsSim, dCompSimSummaryIM, dContSimSummaryIM, dCompSimSummaryMS, dContSimSummaryMS, dCompSimSummaryRMS, dContSimSummaryRMS = summaryResults[simId][simUID]
            _, dCompSimFullIM, dContSimFullIM, dCompSimFullMS, dContSimFullMS, dCompSimFullRMS, dContSimFullRMS = fullResults[simId][simUID]
            
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
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simCompare = SimCompare()
    simCompare.run()
