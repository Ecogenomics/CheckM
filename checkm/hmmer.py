###############################################################################
#
# hmmer.py - runs HMMER and provides functions for parsing output
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
import subprocess
from re import split as re_split


class FormatError(BaseException):
    pass


class HMMERError(BaseException):
    pass


class HMMMERModeError(BaseException):
    pass


class HMMERRunner():
    """Wrapper for running HMMER3."""
    def __init__(self, mode="dom"):
        self.logger = logging.getLogger('timestamp')

        # make sure HMMER is installed
        self.checkForHMMER()

        # set the mode
        if mode == "dom":
            self.mode = 'domtblout'
        elif mode == 'align':
            self.mode = 'align'
        elif mode == 'fetch':
            self.mode = 'fetch'
        else:
            raise HMMMERModeError("Mode %s not understood" % mode)

    def search(self, db, query, tableOut, hmmerOut, cmdlineOptions='', bKeepOutput=True):
        """Run hmmsearch"""
        # make the output dir and files
        if self.mode != 'domtblout':
            raise HMMMERModeError("Mode %s not compatible with search" % self.mode)

        if not bKeepOutput:
            hmmerOut = '/dev/null'

        cmd = ('hmmsearch --%s %s %s %s %s > %s' % (self.mode, tableOut, cmdlineOptions, db, query, hmmerOut))
        rtn = os.system(cmd)
        if rtn != 0:
            self.logger.error(f'hmmsearch exited with code: {rtn}')
            sys.exit(rtn)

    def align(self, db, query, outputFile, writeMode='>', outputFormat='PSIBLAST', trim=True):
        """Run hmmalign"""
        if self.mode != 'align':
            raise HMMMERModeError("Mode %s not compatible with align" % self.mode)

        # build up the option string
        opts = ''
        if trim:
            opts += ' --trim '

        # run hmmer
        cmd = 'hmmalign --allcol %s --outformat %s %s %s %s %s' % (opts, outputFormat, db, query, writeMode, outputFile)
        rtn = os.system(cmd)
        if rtn == 256:
            # assume user has a newer version of HMMER (>= 3.1b1) and the allcol parameter is no longer valid
            cmd = cmd.replace('--allcol', '')
            rtn = os.system(cmd)

        if rtn != 0:
            self.logger.error(f'hmmalign exited with code: {rtn}')
            sys.exit(1)

    def fetch(self, db, key, fetchFileName, bKeyFile=False):
        """Run hmmfetch"""
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)

        keyFileOpt = ''
        if bKeyFile:
            keyFileOpt = '-f'

        rtn = os.system('hmmfetch ' + keyFileOpt + ' %s %s > %s' % (db, key, fetchFileName))
        if rtn != 0:
            self.logger.error(f'hmmfetch exited with code: {rtn}')
            sys.exit(rtn)

    def index(self, hmmModelFile):
        """Index a HMM file."""
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)

        rtn = os.system('hmmfetch --index %s > /dev/null' % hmmModelFile)
        if rtn != 0:
            self.logger.error(f'hmmfetch exited with code: {rtn}')
            sys.exit(rtn)

    def checkForHMMER(self):
        """Check to see if HMMER is on the system path."""
        try:
            subprocess.call(['hmmsearch', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("Make sure HMMER executables (e.g., hmmsearch, hmmfetch) are on your system path.")
            sys.exit(1)


class HMMERParser():
    """Parses tabular output in domtblout format."""
    def __init__(self, fileHandle, mode='dom'):
        self.handle = fileHandle
        if mode == 'dom':
            self.mode = 'domtblout'
        else:
            raise HMMERError("Mode %s not understood, please use 'dom'" % mode)

    def next(self):
        """Process the next hit in the file, returning None if there are no more hits."""
        hit = self.readHitsDOM()
        if hit == {}:
            return None

        return hit

    def readHitsDOM(self):
        """Process single hit in domtblout format."""
        """
We expect the line to look like:
NODE_925902_length_6780_cov_18.428171_754_2 -            399 PGK                  PF00162.14   384  2.2e-164  543.7   0.1   1   1  1.3e-167  2.5e-164  543.5   0.1     1   384     9   386     9   386 1.00 # 1767 # 2963 # -1 # ID=754_2;partial=00;start_type=ATG;rbs_motif=AGGA;rbs_spacer=5-10bp
        """
        while (1):
            line = self.handle.readline().rstrip()
            try:
                if line[0] != '#' and len(line) != 0:
                    dMatch = re_split(r'\s+', line.rstrip())
                    if len(dMatch) < 23:
                        raise FormatError("Error processing line:\n%s" % (line))
                    refined_match = dMatch[0:22] + [" ".join([str(i) for i in dMatch[22:]])]
                    return HmmerHitDOM(refined_match)
            except IndexError:
                return {}


class HmmerHitDOM():
    """Encapsulate a HMMER hit given in domtblout format."""
    def __init__(self, values):
        if len(values) == 23:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.target_length = int(values[2])
            self.query_name = values[3]

            self.query_accession = values[4]
            if self.query_accession == '-':
                self.query_accession = self.query_name

            self.query_length = int(values[5])
            self.full_e_value = float(values[6])
            self.full_score = float(values[7])
            self.full_bias = float(values[8])
            self.dom = int(values[9])
            self.ndom = int(values[10])
            self.c_evalue = float(values[11])
            self.i_evalue = float(values[12])
            self.dom_score = float(values[13])
            self.dom_bias = float(values[14])
            self.hmm_from = int(values[15])
            self.hmm_to = int(values[16])
            self.ali_from = int(values[17])
            self.ali_to = int(values[18])
            self.env_from = int(values[19])
            self.env_to = int(values[20])
            self.acc = float(values[21])
            self.target_description = values[22]

    def __str__(self):
        return "\t".join([self.target_name,
                          self.target_accession,
                          str(self.target_length),
                          self.query_name,
                          self.query_accession,
                          str(self.query_length),
                          str(self.full_e_value),
                          str(self.full_score),
                          str(self.full_bias),
                          str(self.dom),
                          str(self.ndom),
                          str(self.c_evalue),
                          str(self.i_evalue),
                          str(self.dom_score),
                          str(self.dom_bias),
                          str(self.hmm_from),
                          str(self.hmm_to),
                          str(self.ali_from),
                          str(self.ali_to),
                          str(self.env_from),
                          str(self.env_to),
                          str(self.acc),
                          self.target_description]
                         )
