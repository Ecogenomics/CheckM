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
Create file mapping each internal node to the marker set selected by simulation.
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

import dendropy
from  dendropy.dataobject.taxon import Taxon

class CreateSelectedMarkerSetMapping(object):
    def __init__(self):
        pass
    
    def __inferredMarkerSet(self, tree, node, inferredMarkerSet):
        while node != None:
            uniqueId = node.label.split('|')[0] 
            if inferredMarkerSet.get(uniqueId, 'NA') != 'NA':
                return inferredMarkerSet[uniqueId]
            
            # return domain-specific set if reached that far
            if uniqueId == 'UID2':
                return 'UID2'
            elif uniqueId == 'UID203':
                return 'UID203'
            
            node = node.parent_node

        # started at root, so just return its unique id
        print uniqueId
        return uniqueId
        
    def run(self):   
        fout = open('./simulations/selected_marker_sets.tsv', 'w')    
        
        # read tree
        treeFile = os.path.join('/srv/whitlam/bio/db/checkm', 'genome_tree', 'genome_tree_prok.refpkg', 'genome_tree.final.tre')
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', as_rooted=True, preserve_underscores=True)
        
        # reading marker set inferred via simulation
        inferredMarkerSet = {}
        with open('./simulations/simInferBestMarkerSet.tsv') as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                
                uid = lineSplit[0].split('|')[0].strip()
                markerSetId = lineSplit[1].split('|')[0].strip()
                inferredMarkerSet[uid] = markerSetId
                
        # determine selected marker set for all internal nodes
        for node in tree.internal_nodes():
            uniqueId = node.label.split('|')[0] 
            selectedUID = self.__inferredMarkerSet(tree, node, inferredMarkerSet)
            
            fout.write(uniqueId + '\t' + selectedUID + '\n')
                 
        fout.close()
                  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    createSelectedMarkerSetMapping = CreateSelectedMarkerSetMapping()
    createSelectedMarkerSetMapping.run()
