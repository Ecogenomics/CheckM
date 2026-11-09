###############################################################################
#
# placement.py - place genomes in the reference genome tree with gPlace
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
import re
import sys
import stat
import shutil
import subprocess
import logging
from collections import defaultdict

from checkm.defaultValues import DefaultValues

from checkm.common import checkDirExists
from checkm.util.seqUtils import readFasta, writeFasta
from checkm.hmmer import runWithStdout, exitOnFailure


# minimum version of gPlace required by CheckM
MIN_GPLACE_VERSION = (0, 1, 0)


def parseGPlaceVersion(versionOutput):
    """Parse the version from the output of 'gplace --version' (e.g., 'gplace: version 0.1.0 ...').

    Returns a tuple of integers, or None if the version can't be determined.
    """
    m = re.search(r'version\s+(\d+(?:\.\d+)*)', versionOutput)
    if not m:
        return None

    return tuple(int(x) for x in m.group(1).split('.'))


class GPlaceRunner():
    """Wrapper for running gPlace (https://github.com/cmc-aau/gPlace).

    gPlace is a replacement for pplacer that gives pplacer's placements, bit for bit,
    using less memory and time. It takes pplacer's command line ('gplace place') and
    provides an equivalent of 'guppy tog' ('gplace tog').
    """

    def __init__(self, threads):
        self.logger = logging.getLogger('timestamp')
        self.numThreads = threads

        # make sure a suitable version of gPlace is on the system path
        self._checkForGPlace()

    def run(self, binFiles, resultsParser, outDir, bReducedTree):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        treeFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_TREE_OUT)
        placementJsonOut = os.path.join(alignOutputDir, DefaultValues.PPLACER_JSON_OUT)
        placementOut = os.path.join(alignOutputDir, DefaultValues.PPLACER_OUT)

        # create concatenated alignment file for each bin
        concatenatedAlignFile = self._createConcatenatedAlignment(binFiles, resultsParser, alignOutputDir)

        refPkg = DefaultValues.PPLACER_REF_PACKAGE_FULL
        if bReducedTree:
            refPkg = DefaultValues.PPLACER_REF_PACKAGE_REDUCED

        # check if concatenated alignment file is empty
        # (this can occur when all genomes have no phylogenetically informative marker genes)
        if os.stat(concatenatedAlignFile)[stat.ST_SIZE] == 0:
            self.logger.info('No genomes were identified that could be placed in the reference genome tree.')
            shutil.copyfile(os.path.join(refPkg, DefaultValues.GENOME_TREE), treeFile)
            return

        # run gPlace to place bins in reference genome tree
        self.logger.info('Placing %d bins into the genome tree with gPlace.' % len(binFiles))
        cmd = ['gplace', 'place', '-j', str(self.numThreads),
               '-c', refPkg,
               '-o', placementJsonOut,
               concatenatedAlignFile]
        rtn = runWithStdout(cmd, placementOut)
        exitOnFailure(self.logger, 'gplace place', rtn)

        # extract tree
        rtn = subprocess.call(['gplace', 'tog', '-o', treeFile, placementJsonOut])
        exitOnFailure(self.logger, 'gplace tog', rtn)

    def _createConcatenatedAlignment(self, binFiles, resultsParser, alignOutputDir):
        """Create a concatenated alignment of marker genes for each bin."""

        # read alignment files
        self.logger.info('Reading marker alignment files.')
        alignments = defaultdict(dict)
        files = os.listdir(alignOutputDir)
        binIds = set()
        for f in files:
            if f.endswith('.masked.faa'):
                markerId = f[0:f.find('.masked.faa')]
                seqs = readFasta(os.path.join(alignOutputDir, f))

                for seqId, seq in seqs.items():
                    binId = seqId[0:seqId.find(DefaultValues.SEQ_CONCAT_CHAR)]

                    alignments[markerId][binId] = seq
                    binIds.add(binId)

        # get all markers and their lengths
        markerIds = resultsParser.models[list(resultsParser.models.keys())[0]].keys()
        markerIdLens = {}
        for markerId in markerIds:
            markerIdLens[markerId] = resultsParser.models[list(resultsParser.models.keys())[0]][markerId].leng

        # create concatenated alignment
        self.logger.info('Concatenating alignments.')
        concatenatedSeqs = {}
        for markerId in sorted(markerIds):
            seqs = alignments[markerId]

            for binId in sorted(binIds):
                if binId in seqs:
                    # append alignment
                    concatenatedSeqs[binId] = concatenatedSeqs.get(binId, '') + seqs[binId]
                else:
                    # missing gene
                    concatenatedSeqs[binId] = concatenatedSeqs.get(binId, '') + '-' * markerIdLens[markerId]

        # save concatenated alignment
        concatenatedAlignFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_CONCAT_SEQ_OUT)
        writeFasta(concatenatedSeqs, concatenatedAlignFile)

        return concatenatedAlignFile

    def _checkForGPlace(self):
        """Check that gPlace is on the system path and is a supported version."""

        minVersionStr = '.'.join(str(v) for v in MIN_GPLACE_VERSION)
        installStr = ('gPlace >= %s is required; install it with "pip install gtdb-gplace" '
                      '(requires Python >= 3.12; it can be installed in a separate environment '
                      'as long as the gplace command is on your system path).' % minVersionStr)

        try:
            proc = subprocess.run(['gplace', '--version'], stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT)
        except OSError:
            self.logger.error('Make sure gplace is on your system path. ' + installStr)
            sys.exit(1)

        versionOutput = proc.stdout.decode('utf-8', errors='replace').strip()
        version = parseGPlaceVersion(versionOutput)
        if proc.returncode != 0 or version is None:
            self.logger.error("Unable to determine the version of gPlace from 'gplace --version': %s" % versionOutput)
            self.logger.error(installStr)
            sys.exit(1)

        if version < MIN_GPLACE_VERSION:
            self.logger.error('gPlace %s was found, but %s' % ('.'.join(str(v) for v in version), installStr))
            sys.exit(1)

        self.version = version
