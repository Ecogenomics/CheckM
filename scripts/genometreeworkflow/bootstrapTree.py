#!/usr/bin/env python3

"""
Performs bootstrap replicates.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import random

class BootstrapTree(object):
    def __init__(self):
        pass

    def readFasta(self, seqFile):
        try:
            fh = file(seqFile)
        except IOError:
            print(("File '" + seqFile + "' does not exist."))
            sys.exit()

        seqs = {}
        for line in fh:
            if line.startswith('>'):
                seqId = line[1:].rstrip('\n')
                seqs[seqId] = ''
            else:
                seqs[seqId] += line.rstrip('\n').rstrip('*')

        return seqs

    def bootstrap(self, seqs, outputFile):
        alignmentLen = len(seqs[list(seqs.keys())[0]])
        cols = [random.randint(0, alignmentLen-1) for _ in range(alignmentLen)]

        fout = open(outputFile, 'w')
        for seqId in seqs:
            fout.write('>' + seqId + '\n')
            seq = seqs[seqId]
            for col in cols:
                fout.write(seq[col])
            fout.write('\n')

        fout.close()

    def run(self, bootstrapDir, treeFile, seqFile, numBootstraps, numProcessors, outputTree):
        seqs = self.readFasta(seqFile)

        # create bootstrap trees
        print('Creating bootstrap alignments.')
        treeListFile = os.path.join(bootstrapDir, 'bootstrap_trees.txt')
        treeListOut = open(treeListFile, 'w')

        bootstrapTreeFiles = []
        for i in range(0, numBootstraps):
            bootstrapTreeFile = os.path.join(bootstrapDir, 'bootstrap_tree.' + str(i) + '.tre')
            bootstrapTreeFiles.append(bootstrapTreeFile)

            alnFile = os.path.join(bootstrapDir, 'bootstrap_aln.' + str(i) + '.faa')
            self.bootstrap(seqs, alnFile)

            cmd = 'FastTree -quiet -nosupport -wag -gamma ' + alnFile + ' > ' + bootstrapTreeFile + '\n'
            treeListOut.write(cmd)

        treeListOut.close()
        
        print('Building bootstrap trees.')
        os.system('cat ' + treeListFile + ' | parallel --max-procs ' + str(numProcessors))

        # create single file with bootstrap trees
        bootstrapFile = os.path.join(bootstrapDir, 'bootstraps.all.tre')
        bootstrapOut = open(bootstrapFile, 'w')
        for bootstrapTreeFile in bootstrapTreeFiles:
            for line in open(bootstrapTreeFile):
                bootstrapOut.write(line)
        bootstrapOut.close()

        print(('  Bootstrap trees written to: ' + bootstrapFile))

        # determine bootstrap support for original tree
        print('Determining bootstrap support for original tree.')
        os.system('CompareToBootstrap.pl -tree ' + treeFile + ' -boot ' + bootstrapFile + ' > ' + outputTree)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Performs bootstrap replicates.",
                                      epilog="FastTree arguments may also be provided (e.g., -nt, -wag).")

    parser.add_argument('bootstrap_dir', help='directory to store bootstrap alignments and trees')
    parser.add_argument('tree_file', help='original tree requiring bootstrap support values')
    parser.add_argument('seq_file', help='multiple sequence alignment used to build original tree (in FASTA format).')
    parser.add_argument('output_tree', help='output tree with bootstrap values.')
    parser.add_argument('-b', '--num_bootstraps', help='number of bootstrap replicates to perform (default = 100).', type=int, default=100)
    parser.add_argument('-t', '--threads', help='number of threads (default = 1).', type = int, default=1)

    args = parser.parse_args()

    bootstrapTree = BootstrapTree()
    bootstrapTree.run(args.bootstrap_dir, args.tree_file, args.seq_file, args.num_bootstraps, args.threads, args.output_tree)
