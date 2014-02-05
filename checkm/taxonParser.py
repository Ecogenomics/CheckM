###############################################################################
#
# taxonParser.py - parse taxonomic-specific marker sets
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
from collections import defaultdict

import prettytable

from markerSet import MarkerSet
from lib.taxonomyUtils import taxonomicRanks

import defaultValues

class TaxonParser():
    """Parse taxonomic-specific marker sets."""
    def __init__(self):
        self.logger = logging.getLogger()
        
    def __readMarkerSets(self):
        taxonMarkerSetFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'taxon_marker_sets.tsv')
        
        taxonMarkerSets = defaultdict(dict)
        for line in open(taxonMarkerSetFile):
            lineSplit = line.split('\t')
            rank = lineSplit[0]
            taxon = lineSplit[1]
            numGenomes = int(lineSplit[2])
            markerSet = eval(lineSplit[5].rstrip())
            
            taxonMarkerSets[rank][taxon] = [numGenomes, MarkerSet(markerSet)]
            
        return taxonMarkerSets
        
    def list(self, rankFilter='ALL'):
        """ List all available marker sets from the specified rank."""

        taxonMarkerSets = self.__readMarkerSets()
        
        header = ['Rank', 'Taxon', '# genomes', '# marker genes', '# marker sets']
        pTable = prettytable.PrettyTable(header)
        pTable.align = 'c'
        pTable.align['Rank'] = 'l'
        pTable.align['Taxon'] = 'l'
        pTable.hrules = prettytable.FRAME
        pTable.vrules = prettytable.NONE
        
        for rank in taxonomicRanks:
            if rankFilter == 'ALL' or rankFilter == rank:
                for taxon in sorted(taxonMarkerSets[rank]):
                    numGenomes, markerSet = taxonMarkerSets[rank][taxon]
                    
                    numMarkers, numMarkerSets = markerSet.size()
                    pTable.add_row([rank, taxon, numGenomes, numMarkers, numMarkerSets])
               
        print ''     
        print pTable.get_string()
            
    def markerSet(self, rank, taxon, markerFile):
        """Obtain specified taxonomic-specific marker set."""
    
        taxonMarkerSets = self.__readMarkerSets()
        
        if rank not in taxonMarkerSets:
            self.logger.error('  Unrecognized taxonomic rank: ' + rank)
            return False
        elif taxon not in taxonMarkerSets[rank]:
            self.logger.error('  Unrecognized taxon: %s (in rank %s): ' % (taxon, rank))
            return False
        
        numGenomes, markerSet = taxonMarkerSets[rank][taxon]
        numMarkers, numMarkerSets = markerSet.size()
        
        self.logger.info('  Marker set for %s contains %d marker genes arranged in %d sets.' % (taxon, numMarkers, numMarkerSets))
        self.logger.info('    Marker set inferred from %d reference genomes.' % numGenomes)
        
        fout = open(markerFile, 'w')
        fout.write(defaultValues.TAXON_MARKER_FILE_HEADER + '\n')        
        fout.write(rank + '\t' + taxon + '\t' + str(numMarkers) + '\t' + str(numMarkerSets) + '\t' + str(markerSet) + '\n')
        fout.close()
        
        return True
