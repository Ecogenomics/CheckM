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

__prog_desc__ = 'create dictionary of gene counts for each genome'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse

from checkm.lib.img import IMG

class GeneCountDict(object):
    def __init__(self):
        pass

    def run(self, outputFile):
        img = IMG()
        
        print 'Identifying all IMG prokaryotic genomes with valid data.'
        metadata = img.genomeMetadata()
        genomeIds = img.genomeIdsByTaxonomy('prokaryotes', metadata)
        genomeMissingData = img.genomesWithMissingData(genomeIds)
        genomeIds -= genomeMissingData
        
        print '  Identified %d valid genomes.' % (len(genomeIds))
        
        print 'Calculating gene copy number for each genome.'
        countTable = img.geneCountTable(genomeIds)
         
        fout = open(outputFile, 'w')
        fout.write(str(countTable))
        fout.close()
        
        print 'Gene count dictionary to: ' + outputFile
        
if __name__ == '__main__':
    print 'GeneCountDict v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_file', help='output file')

    args = parser.parse_args()

    geneCountDict = GeneCountDict()
    geneCountDict.run(args.output_file)
