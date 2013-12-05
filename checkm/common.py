###############################################################################
#
# common.py - utility functions used in many places in CheckM
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
import errno
import sys
import ast

import numpy as np

def readDistribution(distPer, dirPostfix):
    """Read distribution file."""
    distFile = os.path.join(os.path.dirname(sys.argv[0]) + '/../data/', dirPostfix, 'distribution_' + distPer + '.txt')
    checkFileExists(distFile)
    
    with open(distFile, 'r') as f:
        s = f.read()
        d = ast.literal_eval(s)
        
    return d

def findNearest(array, value):
    '''Find nearest array element to a given value.'''
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def checkFileExists(inputFile):
    if not os.path.exists(inputFile):
        sys.stderr.write('[Error] Input file does not exists: ' + inputFile + '\n')
        sys.exit()

def makeSurePathExists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        
def binIdFromFilename(filename):
    binId = os.path.basename(filename)
    binId = os.path.splitext(binId)[0]
    
    return binId
        
def reassignStdOut(outFile):
    oldStdOut = sys.stdout
    if(outFile != ''):
        try:
            # redirect stdout to a file
            sys.stdout = open(outFile, 'w')
        except:
            print("Error diverting stout to file: ", sys.exc_info()[0])
            raise
            
    return oldStdOut

def restoreStdOut(outFile, oldStdOut, bQuiet=False):
    if(outFile != ''):
        try:
            # redirect stdout to a file
            sys.stdout.close()
            sys.stdout = oldStdOut
        except:
            print("Error restoring stdout ", sys.exc_info()[0])
            raise