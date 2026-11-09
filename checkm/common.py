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
import logging

import numpy as np

from checkm.defaultValues import DefaultValues


def getBinIdsFromOutDir(outDir):
    """Get bin ids."""
    binIds = []

    binDir = os.path.join(outDir, 'bins')

    for f in os.listdir(binDir):
        if os.path.isdir(os.path.join(binDir, f)) and f != 'storage':
            binIds.append(f)

    return binIds


def readDistribution(prefix):
    """Read distribution file."""
    distFile = os.path.join(DefaultValues.DISTRIBUTION_DIR, prefix + '.txt')
    checkFileExists(distFile)

    with open(distFile, 'r') as f:
        s = f.read()
        d = ast.literal_eval(s)

    return d


def findNearest(array, value):
    '''Find nearest array element to a given value.'''

    idx = (np.abs(np.array(array) - value)).argmin()
    return array[idx]


def checkBinInputExists(binInput, bCalledGenes):
    """Check if bin input exists.

    The bin input is either a directory or a tab-separated table. For a
    table, the first non-blank line that does not start with '#' must
    have at least 2 columns (3 columns if bCalledGenes is True). Exits
    with an error if these conditions are not met.
    """
    if os.path.isdir(binInput):
        return

    checkFileExists(binInput)

    first_row = None
    with open(binInput) as binInput_ih:
        for line in binInput_ih:
            if line.strip() and not line.lstrip().startswith('#'):
                first_row = line.strip().split("\t")
                break

    num_required_columns = 3 if bCalledGenes else 2
    if first_row is None or len(first_row) < num_required_columns:
        logger = logging.getLogger("timestamp")
        logger.error('Input table format is not correct: ' +
                     binInput + '\n')
        logger.error(
            'If input is nucleotide contigs, please supply at least 2 columns: [genome_id, genome_contigs_file]')
        logger.error(
            'If input is gene, please supply at least 3 columns: [genome_id, genome_contigs_file, genome_genes_file]')
        sys.exit(1)


def checkFileExists(inputFile):
    """Check if file exists."""
    if not os.path.exists(inputFile):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exist: ' + inputFile + '\n')
        sys.exit(1)


def checkDirExists(inputDir):
    """Check if directory exists."""
    if not os.path.isdir(inputDir):
        logger = logging.getLogger('timestamp')
        logger.error('Input directory does not exist: ' + inputDir + '\n')
        sys.exit(1)


def makeSurePathExists(path):
    """Create directory if it does not exist."""
    if not path:
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path does not exist: ' + path + '\n')
            sys.exit(1)


def binIdFromFilename(filename):
    """Extract bin id from bin filename."""

    binId = os.path.basename(filename)

    # remove common compression extensions
    compress_ext = ['.gz', '.bz2', '.rar', '.tgz', '.zip', '.7z']
    for ext in compress_ext:
        if binId.endswith(ext):
            binId = binId[0:len(binId) - len(ext)]
            break

    binId = os.path.splitext(binId)[0]

    return binId


def reassignStdOut(outFile):
    """Redirect standard out to a file."""
    oldStdOut = sys.stdout
    if(outFile != ''):
        try:
            # redirect stdout to a file
            sys.stdout = open(outFile, 'w')
        except:
            logger = logging.getLogger('timestamp')
            logger.error("Error diverting stdout to file: " + outFile)
            sys.exit(1)

    return oldStdOut


def restoreStdOut(outFile, oldStdOut):
    """Redirect standard out back to system standard out."""
    if(outFile != ''):
        try:
            # restore original stdout
            sys.stdout.close()
            sys.stdout = oldStdOut
        except:
            logger = logging.getLogger('timestamp')
            logger.error("Error restoring stdout: %s", outFile)
            sys.exit(1)
