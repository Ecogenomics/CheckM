###############################################################################
#
# prodigal.py - runs prodigal and provides functions for parsing output
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
import stat
import subprocess
import logging
import shutil
import tempfile

import numpy as np

from checkm.defaultValues import DefaultValues
from checkm.common import checkFileExists
from checkm.util.seqUtils import readFasta, writeFasta


class ProdigalError(BaseException):
    pass


class ProdigalRunner():
    """Wrapper for running prodigal."""

    def __init__(self, outDir):
        self.logger = logging.getLogger('timestamp')

        # make sure prodigal is installed
        self.checkForProdigal()

        self.aaGeneFile = os.path.join(outDir, DefaultValues.PRODIGAL_AA)
        self.ntGeneFile = os.path.join(outDir, DefaultValues.PRODIGAL_NT)
        self.gffFile = os.path.join(outDir, DefaultValues.PRODIGAL_GFF)

    def run(self, query, bNucORFs=True):

        prodigal_input = query

        # gather statistics about query file
        seqs = readFasta(prodigal_input)
        totalBases = 0
        for seqId, seq in seqs.items():
            totalBases += len(seq)

        # decompress gzip input files
        tmp_dir = None
        if prodigal_input.endswith('.gz'):
            tmp_dir = tempfile.mkdtemp()
            prodigal_input = os.path.join(
                tmp_dir, os.path.basename(prodigal_input[0:-3]))
            writeFasta(seqs, prodigal_input)

        # call ORFs with different translation tables and select the one with the highest coding density
        tableCodingDensity = {}
        for translationTable in [4, 11]:
            aaGeneFile = self.aaGeneFile + '.' + str(translationTable)
            ntGeneFile = self.ntGeneFile + '.' + str(translationTable)
            gffFile = self.gffFile + '.' + str(translationTable)

            # check if there is sufficient bases to calculate prodigal parameters
            if totalBases < 100000:
                procedureStr = 'meta'  # use best precalculated parameters
            else:
                procedureStr = 'single'  # estimate parameters from data

            if bNucORFs:
                cmd = ('prodigal -p %s -q -m -f gff -g %d -a %s -d %s -i %s > %s 2> /dev/null' % (procedureStr,
                                                                                                  translationTable,
                                                                                                  aaGeneFile,
                                                                                                  ntGeneFile,
                                                                                                  prodigal_input,
                                                                                                  gffFile))
            else:
                cmd = ('prodigal -p %s -q -m -f gff -g %d -a %s -i %s > %s 2> /dev/null' % (procedureStr,
                                                                                            translationTable,
                                                                                            aaGeneFile,
                                                                                            prodigal_input,
                                                                                            gffFile))


            rtn = os.system(cmd)
            if rtn != 0 or not self._areORFsCalled(aaGeneFile):
                # prodigal will fail to learn a model if the input genome has a large number of N's
                # so try gene prediction with 'meta'
                msg = f"Prodigal failed or returned no output (code: {rtn})."
                if procedureStr == 'single':
                    msg += " Retrying with '-p meta' due to possible high N content in the genome."
                    self.logger.warning(msg)
                    cmd = cmd.replace('-p single', '-p meta')
                    rtn = os.system(cmd)
                    if rtn != 0 or not self._areORFsCalled(aaGeneFile):
                        self.logger.error(f"Prodigal failed again with '-p meta' (code: {rtn}).")
                        sys.exit(rtn)
                else:
                    self.logger.error(msg)
                    sys.exit(rtn)

            # determine coding density
            prodigalParser = ProdigalGeneFeatureParser(gffFile)

            codingBases = 0
            for seqId, seq in seqs.items():
                codingBases += prodigalParser.codingBases(seqId)

            if totalBases != 0:
                codingDensity = float(codingBases) / totalBases
            else:
                codingDensity = 0
            tableCodingDensity[translationTable] = codingDensity

        # determine best translation table
        bestTranslationTable = 11
        if (tableCodingDensity[4] - tableCodingDensity[11] > 0.05) and tableCodingDensity[4] > 0.7:
            bestTranslationTable = 4

        shutil.copyfile(self.aaGeneFile + '.' +
                        str(bestTranslationTable), self.aaGeneFile)
        shutil.copyfile(self.gffFile + '.' +
                        str(bestTranslationTable), self.gffFile)
        if bNucORFs:
            shutil.copyfile(self.ntGeneFile + '.' +
                            str(bestTranslationTable), self.ntGeneFile)

        # clean up redundant prodigal results
        for translationTable in [4, 11]:
            os.remove(self.aaGeneFile + '.' + str(translationTable))
            os.remove(self.gffFile + '.' + str(translationTable))
            if bNucORFs:
                os.remove(self.ntGeneFile + '.' + str(translationTable))

        if tmp_dir:
            shutil.rmtree(tmp_dir)

        return bestTranslationTable

    def _areORFsCalled(self, aaGeneFile):
        return os.path.exists(aaGeneFile) and os.stat(aaGeneFile)[stat.ST_SIZE] != 0

    def areORFsCalled(self, bNucORFs):
        # if requested, check if nucleotide gene sequences have been generated
        if bNucORFs:
            return os.path.exists(self.ntGeneFile) and os.stat(self.ntGeneFile)[stat.ST_SIZE] != 0

        # otherwise, only the amino acid gene sequences are required
        return os.path.exists(self.aaGeneFile) and os.stat(self.aaGeneFile)[stat.ST_SIZE] != 0

    def checkForProdigal(self):
        """Check to see if Prodigal is on the system before we try to run it."""

        # Assume that a successful prodigal -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(
                ['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("Make sure prodigal is on your system path.")
            sys.exit(1)


class ProdigalFastaParser():
    """Parses prodigal FASTA output."""

    def __init__(self):
        pass

    def genePositions(self, filename):
        checkFileExists(filename)

        gp = {}
        for line in open(filename):
            if line[0] == '>':
                lineSplit = line[1:].split()

                geneId = lineSplit[0]
                startPos = int(lineSplit[2])
                endPos = int(lineSplit[4])

                gp[geneId] = [startPos, endPos]

        return gp


class ProdigalGeneFeatureParser():
    """Parses prodigal FASTA output."""

    def __init__(self, filename):
        checkFileExists(filename)

        self.genes = {}
        self.lastCodingBase = {}

        self._parseGFF(filename)

        self.codingBaseMasks = {}
        for seqId in self.genes:
            self.codingBaseMasks[seqId] = self._buildCodingBaseMask(seqId)

    def _parseGFF(self, filename):
        """Parse genes from GFF file."""
        self.translationTable = None
        for line in open(filename):
            if line.startswith('# Model Data') and not self.translationTable:
                lineSplit = line.split(';')
                for token in lineSplit:
                    if 'transl_table' in token:
                        self.translationTable = int(
                            token[token.find('=') + 1:])

            if line[0] == '#' or line.strip() == '"':
                # work around for Prodigal having lines with just a
                # quotation on it when FASTA files have Windows style
                # line endings
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId not in self.genes:
                geneCounter = 0
                self.genes[seqId] = {}
                self.lastCodingBase[seqId] = 0

            geneId = seqId + '_' + str(geneCounter)
            geneCounter += 1

            start = int(lineSplit[3])
            end = int(lineSplit[4])

            self.genes[seqId][geneId] = [start, end]
            self.lastCodingBase[seqId] = max(self.lastCodingBase[seqId], end)

    def _buildCodingBaseMask(self, seqId):
        """Build mask indicating which bases in a sequences are coding."""

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes; indices adjusted
        # to account for GFF file using 1-based indexing
        codingBaseMask = np.zeros(self.lastCodingBase[seqId])
        for pos in self.genes[seqId].values():
            codingBaseMask[pos[0]-1:pos[1]] = 1

        return codingBaseMask

    def codingBases(self, seqId, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end)."""

        # check if sequence has any genes
        if seqId not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end == None:
            end = self.lastCodingBase[seqId]

        return np.sum(self.codingBaseMasks[seqId][start:end])
