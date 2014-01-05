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

__prog_desc__ = 'infer trees from sequence alignments'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import argparse

class MakeTrees(object):
    def __init__(self):
        pass
    
    def run(self, alignDir, outputDir, extension, numThreads):
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
    
        treeList = open(outputDir + '/treeList.txt', 'w')
        files = os.listdir(alignDir)
        for f in files:
            if f.endswith(extension):
                # replace any '*' amino acids with an 'X' as many downstream programs do not like asterisk
                fin = open(os.path.join(alignDir, f))
                data = fin.readlines()
                fin.close()
        
                fout = open(os.path.join(alignDir, f), 'w')
                for line in data:
                    if line[0] != '>':
                        line = line.replace('*', 'X')
                    fout.write(line)
                fout.close()
            
                prefix = f[0:f.rfind('.')]
                cmd = 'FastTree -quiet -nosupport -wag -gamma -log ' + outputDir + '/' + prefix + '.log ' + alignDir + '/' + f + ' > ' + outputDir + '/' + prefix + '.tre'
                treeList.write(cmd + '\n')
    
        treeList.close()
        
        print 'Building trees...'
        os.system('cat ' + outputDir + '/treeList.txt | parallel --max-procs ' + str(numThreads))

if __name__ == '__main__':
    print 'MakeTrees v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'
  
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('align_dir', help='directory containing multiple sequence alignments')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-x', '--extension', help='extension of alignment files to process', default = '.aln.masked.faa')
    parser.add_argument('-t', '--threads', help='number of threads to use', type=int, default = 1)
    
    args = parser.parse_args()
    
    makeTrees = MakeTrees()
    makeTrees.run(args.align_dir, args.output_dir, args.extension, args.threads)
