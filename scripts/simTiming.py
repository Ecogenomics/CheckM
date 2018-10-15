#!/usr/bin/env python3

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
Simulate performance of marker sets under different conditions.
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
import time

import dendropy
from  dendropy.dataobject.taxon import Taxon

from checkm.lib.img import IMG
from lib.markerSetBuilder import MarkerSetBuilder

class Simulation(object):
    def __init__(self):
        self.markerSetBuilder = MarkerSetBuilder()
        self.img = IMG()

    def run(self):
        print('\n  Reading reference genome tree.')
        treeFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # get all Finished, Trusted genomes
        metadata = self.img.genomeMetadata()
        bacteriaIds = self.img.getGenomesByClade('domain', 'Bacteria', metadata)
        print(('# Bacteria: %d' % len(bacteriaIds)))
        
        start = time.time()
        #self.markerSetBuilder.cachedGeneCountTable = self.img.geneCountTable(metadata.keys())
        end = time.time()
        print(('globalGeneCountTable: %.2f' % (end - start)))
        
        start = time.time()
        #self.markerSetBuilder.precomputeGenomeSeqLens(metadata.keys())
        end = time.time()
        print(('precomputeGenomeSeqLens: %.2f' % (end - start)))
        
        start = time.time()
        #self.markerSetBuilder.precomputeGenomeFamilyPositions(metadata.keys(), 5000)
        end = time.time()
        print(('precomputeGenomeFamilyPositions: %.2f' % (end - start)))
        
        #start = time.time()
        #test = self.img.geneDistTable(metadata.keys(), self.markerSetBuilder.globalGeneCountTable.keys(), spacingBetweenContigs=1e6)
        #end = time.time()
        #print 'geneDistTable: %.2f' % (end - start)
        
        #t = raw_input('waiting...')
        
        start = time.time()
        #testGenomeId = archaeaIds.pop()
        #testNode = tree.find_node_with_taxon_label('IMG_' + testGenomeId)
        #binMarkerSets, refinedBinMarkerSet = self.markerSetBuilder.buildBinMarkerSet(tree, testNode.parent_node, 0.97, 0.97, bMarkerSet = True)
        markerSet = self.markerSetBuilder.buildMarkerSet(bacteriaIds, 0.97, 0.97)
        end = time.time()
        print(('buildMarkerSet: %.2f' % (end - start)))
        
 
        print((len(markerSet.markerSet)))
        test = eval("[set(['TIGR01080', 'pfam08071', 'pfam00900']), set(['pfam08069', 'pfam00312']), set(['pfam00276', 'pfam00573', 'pfam00297']), set(['pfam00333', 'pfam00327', 'pfam03719']), set(['pfam04563', 'pfam04560', 'pfam04983', 'pfam04566', 'pfam04567', 'pfam04565', 'pfam04997', 'pfam00562', 'pfam05000', 'pfam00623', 'pfam01191', 'pfam04561', 'TIGR02389']), set(['pfam00831', 'pfam01868', 'pfam00189', 'pfam00237']), set(['pfam01280', 'pfam00861', 'pfam01655']), set(['pfam01194', 'pfam00380', 'pfam00572']), set(['TIGR00425', 'pfam08068']), set(['pfam01157', 'pfam03874']), set(['pfam00416', 'TIGR01018']), set(['pfam00679', 'pfam03764']), set(['TIGR00344', 'TIGR03683']), set(['pfam01193', 'pfam01000']), set(['pfam00750', 'pfam05746']), set(['pfam00935', 'pfam01667']), set(['pfam00867', 'pfam00752']), set(['pfam01172', 'pfam09377']), set(['pfam03950', 'pfam00749']), set(['pfam00181', 'pfam03947']), set(['pfam00687', 'pfam00466']), set(['TIGR03679', 'TIGR00289']), set(['pfam01198', 'pfam01912']), set(['pfam00673', 'pfam00281']), set(['TIGR00134']), set(['pfam00410']), set(['pfam00411']), set(['pfam01090']), set(['pfam01092']), set(['pfam04919']), set(['TIGR00336']), set(['pfam01864']), set(['TIGR00442']), set(['pfam01866']), set(['pfam01780']), set(['TIGR01046']), set(['pfam00318']), set(['pfam00252']), set(['pfam09173']), set(['pfam00238']), set(['pfam01798']), set(['pfam01246']), set(['pfam07541']), set(['pfam00736']), set(['TIGR00522']), set(['pfam01269']), set(['TIGR00329']), set(['pfam01015']), set(['TIGR00392']), set(['pfam00203']), set(['TIGR00398']), set(['pfam01725']), set(['pfam02005']), set(['TIGR00422']), set(['pfam03439']), set(['pfam01351']), set(['pfam01922']), set(['pfam11987']), set(['pfam04127']), set(['TIGR00064']), set(['TIGR00389']), set(['pfam13656']), set(['pfam00298']), set(['TIGR00432']), set(['TIGR03677']), set(['pfam00958']), set(['pfam05221']), set(['pfam00347']), set(['TIGR03685']), set(['pfam03876']), set(['pfam01192']), set(['pfam01984']), set(['pfam00827']), set(['pfam01982']), set(['pfam01981']), set(['TIGR00408']), set(['TIGR00270']), set(['TIGR03665']), set(['pfam02978']), set(['pfam03484']), set(['pfam01201']), set(['TIGR02076']), set(['pfam00832']), set(['pfam00833']), set(['TIGR00419']), set(['pfam00177']), set(['pfam06418']), set(['TIGR00057']), set(['TIGR00549']), set(['pfam13685']), set(['pfam05670']), set(['pfam01849']), set(['TIGR02338']), set(['TIGR00468']), set(['pfam09249']), set(['pfam01287']), set(['pfam00164']), set(['pfam01282']), set(['TIGR03724']), set(['pfam01200']), set(['TIGR02153']), set(['TIGR00670']), set(['pfam00398']), set(['TIGR01213']), set(['pfam06026']), set(['pfam04019']), set(['pfam04010']), set(['pfam00366'])]")
        print((len(test)))
        
        for ms in markerSet.markerSet:
            bMatch = False
            for tms in test:
                if tms == ms:
                    print(ms)
                    print(tms)
                    print('---------')
                    bMatch = True
                    break
                
            if not bMatch:
                print('BOO!')
        
        if str(markerSet.markerSet) == "[set(['TIGR01080', 'pfam08071', 'pfam00900']), set(['pfam08069', 'pfam00312']), set(['pfam00276', 'pfam00573', 'pfam00297']), set(['pfam00333', 'pfam00327', 'pfam03719']), set(['pfam04563', 'pfam04560', 'pfam04983', 'pfam04566', 'pfam04567', 'pfam04565', 'pfam04997', 'pfam00562', 'pfam05000', 'pfam00623', 'pfam01191', 'pfam04561', 'TIGR02389']), set(['pfam00831', 'pfam01868', 'pfam00189', 'pfam00237']), set(['pfam01280', 'pfam00861', 'pfam01655']), set(['pfam01194', 'pfam00380', 'pfam00572']), set(['TIGR00425', 'pfam08068']), set(['pfam01157', 'pfam03874']), set(['pfam00416', 'TIGR01018']), set(['pfam00679', 'pfam03764']), set(['TIGR00344', 'TIGR03683']), set(['pfam01193', 'pfam01000']), set(['pfam00750', 'pfam05746']), set(['pfam00935', 'pfam01667']), set(['pfam00867', 'pfam00752']), set(['pfam01172', 'pfam09377']), set(['pfam03950', 'pfam00749']), set(['pfam00181', 'pfam03947']), set(['pfam00687', 'pfam00466']), set(['TIGR03679', 'TIGR00289']), set(['pfam01198', 'pfam01912']), set(['pfam00673', 'pfam00281']), set(['TIGR00134']), set(['pfam00410']), set(['pfam00411']), set(['pfam01090']), set(['pfam01092']), set(['pfam04919']), set(['TIGR00336']), set(['pfam01864']), set(['TIGR00442']), set(['pfam01866']), set(['pfam01780']), set(['TIGR01046']), set(['pfam00318']), set(['pfam00252']), set(['pfam09173']), set(['pfam00238']), set(['pfam01798']), set(['pfam01246']), set(['pfam07541']), set(['pfam00736']), set(['TIGR00522']), set(['pfam01269']), set(['TIGR00329']), set(['pfam01015']), set(['TIGR00392']), set(['pfam00203']), set(['TIGR00398']), set(['pfam01725']), set(['pfam02005']), set(['TIGR00422']), set(['pfam03439']), set(['pfam01351']), set(['pfam01922']), set(['pfam11987']), set(['pfam04127']), set(['TIGR00064']), set(['TIGR00389']), set(['pfam13656']), set(['pfam00298']), set(['TIGR00432']), set(['TIGR03677']), set(['pfam00958']), set(['pfam05221']), set(['pfam00347']), set(['TIGR03685']), set(['pfam03876']), set(['pfam01192']), set(['pfam01984']), set(['pfam00827']), set(['pfam01982']), set(['pfam01981']), set(['TIGR00408']), set(['TIGR00270']), set(['TIGR03665']), set(['pfam02978']), set(['pfam03484']), set(['pfam01201']), set(['TIGR02076']), set(['pfam00832']), set(['pfam00833']), set(['TIGR00419']), set(['pfam00177']), set(['pfam06418']), set(['TIGR00057']), set(['TIGR00549']), set(['pfam13685']), set(['pfam05670']), set(['pfam01849']), set(['TIGR02338']), set(['TIGR00468']), set(['pfam09249']), set(['pfam01287']), set(['pfam00164']), set(['pfam01282']), set(['TIGR03724']), set(['pfam01200']), set(['TIGR02153']), set(['TIGR00670']), set(['pfam00398']), set(['TIGR01213']), set(['pfam06026']), set(['pfam04019']), set(['pfam04010']), set(['pfam00366'])]":
            print('Good to go!')
        else:
            print('oh, shit!!!!')
      
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    simulation = Simulation()
    simulation.run()
