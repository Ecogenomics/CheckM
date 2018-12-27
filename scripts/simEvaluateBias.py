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
Evaluate bias in completeness and contamination estimates
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

from numpy import mean, std

import dendropy
from  dendropy.dataobject.taxon import Taxon

class SimCompare(object):
    def __init__(self):
        #self.simResultsFile = './simulations/simulation.random_scaffolds.cont_per_by_g2.w_refinement_97.draft.tsv.gz'
        self.simResultsFile = './simulations/simulation.draft.w_refinement_50.reduced_param.tsv.gz'
                

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
        
    def __readFullResults(self, filename):
        results = defaultdict(dict)
        genomeIds = set()
        with gzip.open(filename) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                simId = lineSplit[0] + '-' + lineSplit[1] + '-' + lineSplit[2] + '-' + lineSplit[3]
                genomeIds.add(lineSplit[0])
                uid = lineSplit[5].split('|')[0].strip()
                numDescendants = int(lineSplit[6])
                
                compIM = [float(x) for x in lineSplit[9].rstrip().split(',')]
                contIM = [float(x) for x in lineSplit[10].rstrip().split(',')]
                
                compMS = [float(x) for x in lineSplit[11].rstrip().split(',')]
                contMS = [float(x) for x in lineSplit[12].rstrip().split(',')]
                
                compRMS = [float(x) for x in lineSplit[15].rstrip().split(',')]
                contRMS = [float(x) for x in lineSplit[16].rstrip().split(',')]
                
                if len(lineSplit) == 19:
                    # true comp and cont was recorded
                    trueComp = [float(x) for x in lineSplit[17].rstrip().split(',')]
                    trueCont = [float(x) for x in lineSplit[18].rstrip().split(',')]
                else:
                    # assume true comp and cont is extremely close to target comp and cont
                    trueComp = [float(lineSplit[2])*100 for x in range(len(lineSplit[15].rstrip().split(',')))]
                    trueCont = [float(lineSplit[3])*100 for x in range(len(lineSplit[15].rstrip().split(',')))]               
            
                results[simId][uid] = [numDescendants, compIM, contIM, compMS, contMS, compRMS, contRMS, trueComp, trueCont]
                
        print('    Number of test genomes: ' + str(len(genomeIds)))
        
        return results
    
    def tallyResults(self, results, simId, uid, key, compBias, contBias, trueComps, trueConts, compError, contError, correctCompError, correctContError, correctCompBias, correctContBias):

        _, compIM, contIM, compMS, contMS, compRMS, contRMS, trueComp, trueCont = results[simId][uid]

        compBias[key].extend(compRMS)
        contBias[key].extend(contRMS)
        trueComps[key].extend(trueComp)
        trueConts[key].extend(trueCont)
        
        for i in range(0, len(compRMS)):
            compError[key].append(abs(compRMS[i]))
            contError[key].append(abs(contRMS[i]))
            
            t_comp = trueComp[i]/100.0
            t_cont = trueCont[i]/100.0
            adjTrueComp = t_comp + t_cont * ( 1.0 - t_comp)
            adjTrueCont = t_cont - t_cont * ( 1.0 - t_comp)
            
            estComp = compRMS[i] + trueComp[i]
            estCont = contRMS[i] + trueCont[i]
            correctCompError[key].append(abs(estComp - adjTrueComp*100))
            correctContError[key].append(abs(estCont - adjTrueCont*100))
            
            correctCompBias[key].append(estComp - adjTrueComp*100)
            correctContBias[key].append(estCont - adjTrueCont*100)
            

    def run(self):
        print('\n  Reading reference genome tree.')
        treeFile = os.path.join('/srv/whitlam/bio/db/checkm/genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # read simulation results        
        print('  Reading full simulation results.')
        results = self.__readFullResults(self.simResultsFile)
                
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
        print('  Evaluating results.')
        
        compBias = defaultdict(list)
        contBias = defaultdict(list)
        trueComps = defaultdict(list)
        trueConts = defaultdict(list)
        
        compError = defaultdict(list)
        contError = defaultdict(list)
        correctCompError = defaultdict(list)
        correctContError = defaultdict(list)
        
        correctCompBias = defaultdict(list)
        correctContBias = defaultdict(list)
        
        compBiasLineage = defaultdict(list)
        contBiasLineage  = defaultdict(list)
        trueCompsLineage  = defaultdict(list)
        trueContsLineage  = defaultdict(list)
        
        compErrorLineage  = defaultdict(list)
        contErrorLineage  = defaultdict(list)
        correctCompErrorLineage  = defaultdict(list)
        correctContErrorLineage  = defaultdict(list)
        
        correctCompBiasLineage  = defaultdict(list)
        correctContBiasLineage  = defaultdict(list)
        
        
        
        itemsProcessed = 0
        for simId in results:
            itemsProcessed += 1
            statusStr = '  Finished processing %d of %d (%.2f%%) test cases.' % (itemsProcessed, len(results), float(itemsProcessed)*100/len(results))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            genomeId, _, comp, cont = simId.split('-')
            
            key = comp + '-' + cont
                          
            # get results for domain-level marker set
            if 'UID2' in results[simId]:
                uid = 'UID2'
            else:
                uid = 'UID203'
                
            self.tallyResults(results, simId, uid, key, 
                                compBias, contBias, 
                                trueComps, trueConts, 
                                compError, contError, 
                                correctCompError, correctContError, 
                                correctCompBias, correctContBias)
                
            # get results for selected marker set
            simUID = self.__inferredMarkerSet(tree, genomeId, inferredMarkerSet)
            self.tallyResults(results, simId, simUID, key, compBiasLineage, contBiasLineage, 
                                trueCompsLineage, trueContsLineage, 
                                compErrorLineage, contErrorLineage, 
                                correctCompErrorLineage, correctContErrorLineage, 
                                correctCompBiasLineage, correctContBiasLineage)
                            
        sys.stdout.write('\n')
        sys.stdout.flush()
        
        for key in trueComps:
            comp, cont = key.split('-')
            
            print('')
            print('Target comp: %s' % comp)
            print('Target cont: %s' % cont)
            
            meanTrueComps = mean(trueComps[key])/100.0
            meanTrueConts = mean(trueConts[key])/100.0
            print('')
            print('Mean actual comp: %.2f' % (meanTrueComps * 100))
            print('Mean actual cont: %.2f' % (meanTrueConts * 100))
            
            estComp = meanTrueComps + meanTrueConts * ( 1.0 - meanTrueComps)
            estCont = meanTrueConts - meanTrueConts * ( 1.0 - meanTrueComps)
            print('\n------')
            print('Predicted comp bias: %.2f' % ((estComp - meanTrueComps) * 100))
            print('Predicted cont bias: %.2f' % ((estCont - meanTrueConts) * 100))
                
            print('\n------')
            print('Estimated comp bias: %.2f +/- %.2f' % (mean(compBias[key]), std(compBias[key])))
            print('Estimated cont bias: %.2f +/- %.2f' % (mean(contBias[key]), std(contBias[key])))
            print('Estimated comp bias lineage: %.2f +/- %.2f' % (mean(compBiasLineage[key]), std(compBiasLineage[key])))
            print('Estimated cont bias lineage: %.2f +/- %.2f' % (mean(contBiasLineage[key]), std(contBiasLineage[key])))
            
            print('\n------')
            print('Corrected comp bias: %.2f +/- %.2f' % (mean(correctCompBias[key]), std(correctCompBias[key])))
            print('Corrected cont bias: %.2f +/- %.2f' % (mean(correctContBias[key]), std(correctContBias[key])))
            print('Corrected comp bias lineage: %.2f +/- %.2f' % (mean(correctCompBiasLineage[key]), std(correctCompBiasLineage[key])))
            print('Corrected cont bias lineage: %.2f +/- %.2f' % (mean(correctContBiasLineage[key]), std(correctContBiasLineage[key])))
            
            print('\n------')
            print('Abs. error comp: %.2f +/- %.2f' % (mean(compError[key]), std(compError[key])))
            print('Abs. error cont: %.2f +/- %.2f' % (mean(contError[key]), std(contError[key])))
            print('Abs. error comp lineage: %.2f +/- %.2f' % (mean(compErrorLineage[key]), std(compErrorLineage[key])))
            print('Abs. error cont lineage: %.2f +/- %.2f' % (mean(contErrorLineage[key]), std(contErrorLineage[key])))
            
            print('\n------')
            print('Corrected abs. error comp: %.2f +/- %.2f' % (mean(correctCompError[key]), std(correctCompError[key])))
            print('Corrected abs. error cont: %.2f +/- %.2f' % (mean(correctContError[key]), std(correctContError[key])))
            print('Corrected abs. error comp lineage: %.2f +/- %.2f' % (mean(correctCompErrorLineage[key]), std(correctCompErrorLineage[key])))
            print('Corrected abs. error cont lineage: %.2f +/- %.2f' % (mean(correctContErrorLineage[key]), std(correctContErrorLineage[key])))
            
            print('')
            print('*******************************************************************')
            
            fout = open('./simulations/correct_comp_cont_' + key + '.test.tsv','w')
            fout.write('\t'.join(map(str, compBias[key])) + '\n')
            fout.write('\t'.join(map(str, contBias[key])) + '\n')
            fout.write('\t'.join(map(str, compBiasLineage[key])) + '\n')
            fout.write('\t'.join(map(str, contBiasLineage[key])) + '\n')
            
            fout.write('\t'.join(map(str, correctCompBias[key])) + '\n')
            fout.write('\t'.join(map(str, correctContBias[key])) + '\n')
            fout.write('\t'.join(map(str, correctCompBiasLineage[key])) + '\n')
            fout.write('\t'.join(map(str, correctContBiasLineage[key])) + '\n')
            fout.close()

                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simCompare = SimCompare()
    simCompare.run()
