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
import subprocess
import logging
from collections import defaultdict

from checkm.defaultValues import DefaultValues

from checkm.common import checkDirExists
from checkm.util.seqUtils import readFasta, writeFasta
from checkm.hmmer import runWithStdout, exitOnFailure, returnCodeStr
from checkm.util.compression import compressFile, hasContent, uncompressedName


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


def parseMemorySize(value):
    """Convert a memory size for gPlace's --max_memory option to bytes.

    Accepts the same values as gPlace: e.g., '64G', '64GB', '64GiB', '500M', '1.5T',
    or a plain number (GiB). Units are binary (K = 1024, M = 1024^2, ...).
    Raises ValueError if the value is not a valid memory size.
    """
    m = re.fullmatch(r'\s*([0-9]*\.?[0-9]+)\s*([KMGT]?)(I?B)?\s*', str(value), re.IGNORECASE)
    if not m:
        raise ValueError('invalid memory size "%s" (examples: 64G, 500M, 1.5T)' % value)

    unit = m.group(2).upper() or ('' if m.group(3) else 'G')
    size = float(m.group(1)) * 1024 ** ('KMGT'.index(unit) + 1 if unit else 0)
    if size < 1:
        raise ValueError('memory size must be positive: "%s"' % value)

    return int(size)


def gplaceErrors(logFile):
    """Error messages in a gPlace log (e.g., the standard output of 'gplace place')."""
    errors = []
    try:
        with open(logFile) as f:
            for line in f:
                m = re.search(r'\bERROR:\s*(.*)', line)
                if m and not m.group(1).startswith('Controlled exit'):
                    errors.append(m.group(1).strip())
    except OSError:
        pass

    return errors


class GPlaceRunner():
    """Wrapper for running gPlace (https://github.com/cmc-aau/gPlace).

    gPlace is a replacement for pplacer that gives pplacer's placements, bit for bit,
    using less memory and time. It takes pplacer's command line ('gplace place') and
    provides an equivalent of 'guppy tog' ('gplace tog').
    """

    def __init__(self, threads, maxMemory=None):
        """Initialization.

        Parameters
        ----------
        threads : int
            Number of threads used by gPlace.
        maxMemory : str
            Memory gPlace may use (e.g., '64G'), passed to 'gplace place --max_memory'.
            If the placement needs more, gPlace places bins tile by tile within this
            limit, which gives identical results. None leaves the choice to gPlace.
        """
        self.logger = logging.getLogger('timestamp')
        self.numThreads = threads

        self.maxMemory = maxMemory
        if self.maxMemory is not None:
            parseMemorySize(self.maxMemory)  # raises ValueError for invalid sizes

        # make sure a suitable version of gPlace is on the system path
        self._checkForGPlace()

    def run(self, binFiles, resultsParser, outDir):
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

        # check if concatenated alignment file is empty
        # (this can occur when all genomes have no phylogenetically informative marker genes)
        if not hasContent(concatenatedAlignFile):
            self.logger.info('No genomes were identified that could be placed in the reference genome tree.')
            compressFile(os.path.join(refPkg, DefaultValues.GENOME_TREE), treeFile, remove=False)
            return

        # gPlace reads the gzip-compressed alignment, but writes uncompressed
        # files, which are compressed once the tree has been created
        placementJsonOut, compressedPlacementJsonOut = uncompressedName(placementJsonOut), placementJsonOut
        treeFile, compressedTreeFile = uncompressedName(treeFile), treeFile

        # run gPlace to place bins in reference genome tree
        if self.maxMemory is not None:
            self.logger.info('Placing %d bins into the genome tree with gPlace (maximum memory: %s).' % (len(binFiles), self.maxMemory))
        else:
            self.logger.info('Placing %d bins into the genome tree with gPlace.' % len(binFiles))
        cmd = self._placeCmd(refPkg, placementJsonOut, concatenatedAlignFile)
        rtn = runWithStdout(cmd, placementOut)
        if rtn != 0:
            # report gPlace's error message (e.g., --max_memory is too small)
            msg = 'gplace place failed with %s' % returnCodeStr(rtn)
            errors = gplaceErrors(placementOut)
            if errors:
                msg += ': ' + ' '.join(errors)
            else:
                msg += '.'
            self.logger.error(msg)
            self.logger.error('See %s for details.' % placementOut)
            sys.exit(1)

        # extract tree
        rtn = subprocess.call(['gplace', 'tog', '-o', treeFile, placementJsonOut])
        exitOnFailure(self.logger, 'gplace tog', rtn)

        compressFile(placementJsonOut, compressedPlacementJsonOut)
        compressFile(treeFile, compressedTreeFile)

    def _placeCmd(self, refPkg, placementJsonOut, concatenatedAlignFile):
        """Command for placing the concatenated alignment with 'gplace place'."""
        cmd = ['gplace', 'place', '-j', str(self.numThreads),
               '-c', refPkg,
               '-o', placementJsonOut]
        if self.maxMemory is not None:
            cmd += ['--max_memory', self.maxMemory]
        cmd.append(concatenatedAlignFile)

        return cmd

    def _createConcatenatedAlignment(self, binFiles, resultsParser, alignOutputDir):
        """Create a concatenated alignment of marker genes for each bin."""

        # read alignment files (gzip-compressed, or uncompressed if written by CheckM <= 1.2.x)
        self.logger.info('Reading marker alignment files.')
        alignments = defaultdict(dict)
        files = os.listdir(alignOutputDir)
        binIds = set()
        for f in files:
            if f.endswith('.masked.faa') and f + '.gz' in files:
                continue  # use the compressed alignment written by this version of CheckM

            if f.endswith('.masked.faa') or f.endswith('.masked.faa.gz'):
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
