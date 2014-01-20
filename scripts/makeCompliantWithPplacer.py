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

__prog_desc__ = "reformats a tree output by tax2tree's level to make it compliant with pplacer"

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse
import dendropy

class MakeCompliantWithPplacer(object):
    def __init__(self):
        pass

    def run(self, inputTree, outputTree):
        tree = dendropy.Tree.get_from_path(inputTree, schema='newick', as_rooted=True)
        
        for node in tree.internal_nodes():  
            if node.label and ':' in node.label: 
                node.label = 'bs|' + node.label.replace(':', '|').replace(' ', '')
                
        tree.write_to_path(outputTree, schema='newick', suppress_rooting=True)
        
if __name__ == '__main__':
    print 'MakeComplianWithPplacer v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_tree', help='input tree')
    parser.add_argument('output_tree', help='output tree')

    args = parser.parse_args()

    makeCompliantWithPplacer = MakeCompliantWithPplacer()
    makeCompliantWithPplacer.run(args.input_tree, args.output_tree)
