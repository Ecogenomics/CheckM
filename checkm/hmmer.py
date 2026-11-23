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

import io
import os
import sys
import signal
import logging
import contextlib
import subprocess

import pyhmmer

from checkm.util.compression import openFile, compressFile, uncompressedName


class FormatError(BaseException):
    pass


class HMMERError(BaseException):
    pass


class HMMMERModeError(BaseException):
    pass


def returnCodeStr(returnCode):
    """Describe the return code of a subprocess (e.g., 'exit code 1' or 'signal SIGKILL').

    Return codes are those reported by the subprocess module, where a
    negative value -N indicates the program was terminated by signal N.
    """
    if returnCode < 0:
        try:
            return f'signal {signal.Signals(-returnCode).name}'
        except ValueError:
            return f'signal {-returnCode}'

    return f'exit code {returnCode}'


def exitOnFailure(logger, program, returnCode):
    """Log an error and exit with status 1 if an external program failed."""
    if returnCode != 0:
        logger.error(f'{program} failed with {returnCodeStr(returnCode)}.')
        sys.exit(1)


def runWithStdout(cmd, outputFile, append=False):
    """Run a command, given as a list of arguments, writing its standard output to a file.

    Arguments are passed directly to the program without a shell, so
    file paths do not need to be quoted. Returns the subprocess return code.
    """
    with open(outputFile, 'a' if append else 'w') as fout:
        return subprocess.call(cmd, stdout=fout)


def _exitOnHmmerError(logger, program, error):
    """Log an error raised by pyhmmer and exit with status 1."""
    logger.error(f'{program} (pyhmmer) failed: {error!r}')
    sys.exit(1)


def _readProteins(seqFile):
    """Read a protein FASTA file (which may be gzip-compressed) into memory in digital form."""
    with openFile(seqFile, 'rb') as fin:
        data = fin.read()

    with pyhmmer.easel.SequenceFile(io.BytesIO(data), format='fasta', digital=True,
                                    alphabet=pyhmmer.easel.Alphabet.amino()) as seqFileHandle:
        return seqFileHandle.read_block()


def _tableColumn(hits, tableFormat, column):
    """Values of a column of the tabular output of HMMER, in the order hits are reported."""
    table = io.BytesIO()
    hits.write(table, format=tableFormat, header=False)
    return [line.split()[column] for line in table.getvalue().decode().splitlines() if line and line[0] != '#']


def _writeSearchReport(fout, hits):
    """Write the hits to a query HMM as in the main output of hmmsearch, including alignments.

    This follows the layout of hmmsearch (with --notextw). The expected
    number of domains ('exp') and the accuracy of each domain ('acc') are
    not provided by pyhmmer objects, so they are taken from its tabular
    output, which reports hits and domains in the same order. The filter
    statistics and run times of the pipeline are not available.
    """
    hmm = hits.query
    reported = [hit for hit in hits if hit.reported]
    expDomains = _tableColumn(hits, 'targets', 10)
    domainAcc = _tableColumn(hits, 'domains', 21)
    if len(expDomains) != len(reported):
        expDomains = ['-'] * len(reported)
    if len(domainAcc) != sum(1 for hit in reported for d in hit.domains if d.reported):
        domainAcc = None

    def text(value):
        return value.decode() if isinstance(value, bytes) else value

    fout.write('Query:       %s  [M=%d]\n' % (text(hmm.name), hmm.M))
    if hmm.accession:
        fout.write('Accession:   %s\n' % text(hmm.accession))
    if hmm.description:
        fout.write('Description: %s\n' % text(hmm.description))

    # scores for complete sequences
    nameWidth = max([8] + [len(text(hit.name)) for hit in hits])  # as hmmsearch, hits not reported are included
    fout.write('Scores for complete sequences (score includes all domains):\n')
    fout.write('  %22s  %22s  %8s\n' % (' --- full sequence ---', ' --- best 1 domain ---', '-#dom-'))
    fout.write('  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n' % ('E-value', ' score', ' bias', 'E-value', ' score', ' bias',
                                                                    '  exp', 'N', nameWidth, 'Sequence', 'Description'))
    fout.write('  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n' % ('-------', '------', '-----', '-------', '------', '-----',
                                                                    ' ----', '--', nameWidth, '--------', '-----------'))
    bIncludedThreshold = False
    for hit, exp in zip(reported, expDomains):
        if not hit.included and not bIncludedThreshold:
            fout.write('  ------ inclusion threshold ------\n')
            bIncludedThreshold = True

        best = hit.best_domain
        fout.write('  %9.2g %6.1f %5.1f  %9.2g %6.1f %5.1f  %5s %2d  %-*s  %s\n' % (
            hit.evalue, hit.score, hit.bias, best.i_evalue, best.score, best.bias, exp,
            sum(1 for d in hit.domains if d.reported), nameWidth, text(hit.name), text(hit.description) or ''))

    if not reported:
        fout.write('\n   [No hits detected that satisfy reporting thresholds]\n')

    # domains and alignments of each sequence
    fout.write('\n\nDomain annotation for each sequence (and alignments):\n')
    domainIndex = 0
    for hit in reported:
        fout.write('>> %s  %s\n' % (text(hit.name), text(hit.description) or ''))
        domains = [d for d in hit.domains if d.reported]
        if not domains:
            fout.write('   [No individual domains that satisfy reporting thresholds (although complete target did)]\n\n')
            continue

        fout.write(' %3s   %6s %5s %9s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n' % (
            '#', 'score', 'bias', 'c-Evalue', 'i-Evalue', 'hmmfrom', 'hmm to', '  ',
            'alifrom', 'ali to', '  ', 'envfrom', 'env to', '  ', 'acc'))
        fout.write(' %3s   %6s %5s %9s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n' % (
            '---', '------', '-----', '---------', '---------', '-------', '-------', '  ',
            '-------', '-------', '  ', '-------', '-------', '  ', '----'))
        for num, d in enumerate(domains, 1):
            ali = d.alignment
            fout.write(' %3d %s %6.1f %5.1f %9.2g %9.2g %7d %7d %s%s %7d %7d %s%s %7d %7d %s%s %4s\n' % (
                num, '!' if d.included else '?', d.score, d.bias, d.c_evalue, d.i_evalue,
                ali.hmm_from, ali.hmm_to, '[' if ali.hmm_from == 1 else '.', ']' if ali.hmm_to == hmm.M else '.',
                ali.target_from, ali.target_to, '[' if ali.target_from == 1 else '.', ']' if ali.target_to == hit.length else '.',
                d.env_from, d.env_to, '[' if d.env_from == 1 else '.', ']' if d.env_to == hit.length else '.',
                domainAcc[domainIndex + num - 1] if domainAcc else '-'))

        fout.write('\n  Alignments for each domain:\n')
        for num, d in enumerate(domains, 1):
            fout.write('  == domain %d  score: %.1f bits;  conditional E-value: %.2g\n' % (num, d.score, d.c_evalue))
            fout.write(str(d.alignment).rstrip('\n') + '\n\n')
        domainIndex += len(domains)

    if not reported:
        fout.write('\n   [No targets detected that satisfy reporting thresholds]\n')

    fout.write('\n\nInternal pipeline statistics summary:\n')
    fout.write('-------------------------------------\n')
    fout.write('Query model(s):                   %10d  (%d nodes)\n' % (hits.searched_models, hits.searched_nodes))
    fout.write('Target sequences:                 %10d  (%d residues searched)\n' % (hits.searched_sequences, hits.searched_residues))
    fout.write('Initial search space (Z):         %10d\n' % hits.Z)
    fout.write('Domain search space  (domZ):      %10d\n' % hits.domZ)
    fout.write('//\n')


class HMMERRunner():
    """Wrapper for running HMMER3 through pyhmmer.

    pyhmmer provides the HMMER library, so the HMMER programs do not need
    to be installed. The methods give the same results as hmmsearch,
    hmmalign and hmmfetch, which were run by CheckM <= 1.3.2.
    """
    def __init__(self, mode="dom"):
        self.logger = logging.getLogger('timestamp')

        # set the mode
        if mode == "dom":
            self.mode = 'domtblout'
        elif mode == 'align':
            self.mode = 'align'
        elif mode == 'fetch':
            self.mode = 'fetch'
        else:
            raise HMMMERModeError("Mode %s not understood" % mode)

    def search(self, db, query, tableOut, hmmerOut, evalue=0.1, domEvalue=0.1, bKeepOutput=True, threads=1):
        """Search protein sequences with each HMM in a file, as done by hmmsearch.

        The table of domain hits (hmmsearch --domtblout) is written to
        tableOut. If bKeepOutput is True, a report with the alignment of
        each hit (see _writeSearchReport) is written to hmmerOut. Output
        files with a '.gz' extension are gzip-compressed once the search
        has finished. The query may be gzip-compressed. The sequences are
        read once and searched with each HMM, whereas hmmsearch reads the
        sequence file again for each HMM.
        """
        if self.mode != 'domtblout':
            raise HMMMERModeError("Mode %s not compatible with search" % self.mode)

        compressedOutputs = []
        if tableOut.endswith('.gz'):
            compressedOutputs.append((uncompressedName(tableOut) + '.tmp', tableOut))
            tableOut = compressedOutputs[-1][0]
        if bKeepOutput and hmmerOut.endswith('.gz'):
            compressedOutputs.append((uncompressedName(hmmerOut) + '.tmp', hmmerOut))
            hmmerOut = compressedOutputs[-1][0]

        try:
            seqs = _readProteins(query)
            with contextlib.ExitStack() as stack:
                hmmFile = stack.enter_context(pyhmmer.plan7.HMMFile(db))
                tableHandle = stack.enter_context(open(tableOut, 'wb'))
                reportHandle = stack.enter_context(open(hmmerOut, 'w')) if bKeepOutput else None

                if reportHandle:
                    reportHandle.write('\n')

                bHeader = True
                for hits in pyhmmer.hmmsearch(hmmFile, seqs, cpus=threads, E=evalue, domE=domEvalue):
                    hits.write(tableHandle, format='domains', header=bHeader)
                    bHeader = False
                    if reportHandle:
                        _writeSearchReport(reportHandle, hits)

                if reportHandle:
                    reportHandle.write('[ok]\n')
        except Exception as e:
            _exitOnHmmerError(self.logger, 'hmmsearch', e)

        for tmpFile, outputFile in compressedOutputs:
            compressFile(tmpFile, outputFile)

    def align(self, db, query, outputFile, writeMode='>', outputFormat='PSIBLAST', trim=True):
        """Align sequences to the (first) HMM in a file, as done by hmmalign."""
        if self.mode != 'align':
            raise HMMMERModeError("Mode %s not compatible with align" % self.mode)

        try:
            with pyhmmer.plan7.HMMFile(db) as hmmFile:
                hmm = hmmFile.read()
            if hmm is None:
                raise ValueError(f'no HMM in {db}')

            # markers are aligned in parallel by worker processes, so each alignment uses a single thread
            msa = pyhmmer.hmmalign(hmm, _readProteins(query), cpus=1, trim=trim)
            with open(outputFile, 'ab' if writeMode == '>>' else 'wb') as fout:
                msa.write(fout, outputFormat.lower())
        except Exception as e:
            _exitOnHmmerError(self.logger, 'hmmalign', e)

    def fetch(self, db, key, fetchFileName, bKeyFile=False):
        """Write the HMM(s) with the given name or accession to a file, as done by hmmfetch.

        If bKeyFile is True, key is a file with one name or accession per
        line and HMMs are written in this order. HMMs are located with the
        SSI index of the HMM file (<db>.ssi) if it exists, and by reading
        the HMM file otherwise.
        """
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)

        try:
            if bKeyFile:
                with open(key) as f:
                    keys = [line.split()[0] for line in f if line.strip() and line[0] != '#']
            else:
                keys = [key]

            with open(fetchFileName, 'wb') as fout:
                if os.path.exists(db + '.ssi'):
                    self._fetchIndexed(db, keys, fout)
                else:
                    self._fetchUnindexed(db, keys, fout)
        except Exception as e:
            _exitOnHmmerError(self.logger, 'hmmfetch', e)

    def _fetchIndexed(self, db, keys, fout):
        """Write HMMs located with the SSI index of a HMM file."""
        with pyhmmer.easel.SSIReader(db + '.ssi') as index, open(db, 'rb') as fin:
            for key in keys:
                try:
                    fin.seek(index.find_name(key).record_offset)
                except KeyError:
                    raise KeyError(f'HMM {key} not found in {db}') from None

                record = []
                for line in fin:
                    record.append(line)
                    if line.rstrip() == b'//':
                        break

                with pyhmmer.plan7.HMMFile(io.BytesIO(b''.join(record))) as hmmFile:
                    hmmFile.read().write(fout)

    def _fetchUnindexed(self, db, keys, fout):
        """Write HMMs found by reading a HMM file (the first HMM with a given name or accession)."""
        keyToHmm = dict.fromkeys(keys)
        with pyhmmer.plan7.HMMFile(db) as hmmFile:
            for hmm in hmmFile:
                for key in (hmm.name, hmm.accession):
                    key = key.decode() if isinstance(key, bytes) else key
                    if key in keyToHmm and keyToHmm[key] is None:
                        keyToHmm[key] = hmm

        for key in keys:
            if keyToHmm[key] is None:
                raise KeyError(f'HMM {key} not found in {db}')
            keyToHmm[key].write(fout)

    def index(self, hmmModelFile):
        """Create the SSI index of a HMM file (<file>.ssi), as done by 'hmmfetch --index'."""
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)

        try:
            with pyhmmer.easel.SSIWriter(hmmModelFile + '.ssi', exclusive=False) as index, open(hmmModelFile, 'rb') as fin:
                fd = index.add_file(os.path.basename(hmmModelFile), format=0)
                offset = recordOffset = 0
                name = None
                accessions = []
                for line in fin:
                    if line.startswith(b'NAME '):
                        name = line.split()[1].decode()
                    elif line.startswith(b'ACC '):
                        accessions.append(line.split()[1].decode())
                    elif line.rstrip() == b'//':
                        index.add_key(name, fd, recordOffset, 0, 0)
                        for accession in accessions:
                            index.add_alias(accession, name)
                        recordOffset = offset + len(line)
                        name = None
                        accessions = []
                    offset += len(line)
        except Exception as e:
            _exitOnHmmerError(self.logger, 'hmmfetch --index', e)

    def checkForHMMER(self):
        """HMMER is provided by pyhmmer, which is imported with this module."""
        pass


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
                    dMatch = line.split()
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
