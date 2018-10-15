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
Gather all marker gene HMMs into a single model file.
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
import argparse
import tempfile
import uuid

from checkm.taxonParser import TaxonParser
from checkm.treeParser import TreeParser
from checkm.markerSets import MarkerSet
from checkm.hmmer import HMMERRunner
from checkm.lib.pfam import PFAM

class GetMarkerGeneHMMs(object):
    def __init__(self):
        self.outputHMMs = os.path.join('..', 'data', 'hmms', 'checkm.hmm')
        self.hmms = os.path.join('/srv/db/checkm/', 'pfam_tigrfam.hmm')

    def run(self):
        # read all taxonomic-specific marker genes
        print('Reading taxonomic-specific marker genes.')
        taxonomicMarkers = set()
        taxonParser = TaxonParser()
        taxonMarkerSets = taxonParser.readMarkerSets()
        for _, taxa in list(taxonMarkerSets.items()):
            for _, markerSet in list(taxa.items()):
                taxonomicMarkers = taxonomicMarkers.union(markerSet.getMarkerGenes())
                
        print(('  Taxonomic-specific marker genes: %d' % len(taxonomicMarkers)))
                
        # read all lineage-specific marker genes
        print('Reading lineage-specific marker genes.')
        lineageMarkers = set()
        treeParser = TreeParser()
        uniqueIdToLineageStatistics = treeParser.readNodeMetadata()
        for uniqueId, d in list(uniqueIdToLineageStatistics.items()):
            markerSet = MarkerSet(uniqueId, 'NA', int(d['# genomes']), eval(d['marker set']))
            lineageMarkers = lineageMarkers.union(markerSet.getMarkerGenes())
            
        print(('  Lineage-specific marker genes: %d' % len(lineageMarkers)))
        
        # gather all marker genes
        markerGenes = taxonomicMarkers.union(lineageMarkers)
        print(('  Total marker genes: %d' % len(markerGenes)))
        
        # get genes from same clan as marker genes
        print('Gathering HMMs from the same clan as marker genes.')
        pfam = PFAM()
        genesInSameClan = pfam.genesInSameClan(markerGenes)
        allMarkers = markerGenes.union(genesInSameClan)
        
        # create file with all model accession numbers
        keyFile = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        fout = open(keyFile, 'w')
        for modelAcc in allMarkers:
            fout.write(modelAcc + '\n')
        fout.close()
        
        # fetch specified models
        HF = HMMERRunner(mode='fetch')
        HF.fetch(self.hmms, keyFile, self.outputHMMs, bKeyFile=True)
        
        # index the HMM file
        if os.path.exists(self.outputHMMs + '.ssi'):
            os.remove(self.outputHMMs + '.ssi')
        HF.index(self.outputHMMs)
        
        # remove key file
        os.remove(keyFile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Gather all marker gene HMMs into a single model file.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    args = parser.parse_args()

    getMarkerGeneHMMs = GetMarkerGeneHMMs()
    getMarkerGeneHMMs.run()
