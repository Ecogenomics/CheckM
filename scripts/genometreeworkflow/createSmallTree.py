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

__prog_desc__ = 'create reduced tree in order to minimize memory requirements'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import random
import itertools

from checkm.util.seqUtils import readFasta, writeFasta
from checkm.util.img import IMG

import dendropy

class CreateSamllTree(object):
    def __init__(self, outputDir):
        self.__checkForFastTree()

        self.derepConcatenatedAlignFile = os.path.join(outputDir, 'genome_tree.concatenated.derep.fasta')
        self.tree = os.path.join(outputDir, 'genome_tree.final.tre')

        self.img = IMG('/srv/whitlam/bio/db/checkm/img/img_metadata.tsv', '/srv/whitlam/bio/db/checkm/pfam/tigrfam2pfam.tsv')
        self.metadata = self.img.genomeMetadata()

    def __checkForFastTree(self):
        """Check to see if FastTree is on the system path."""

        try:
            exit_status = os.system('FastTree 2> /dev/null')
        except:
            print("Unexpected error!", sys.exc_info()[0])
            raise

        if exit_status != 0:
            print("[Error] FastTree is not on the system path")
            sys.exit()

    def __nearlyIdentical(self, string1, string2, max_diff_perc=0.08):
        max_diff = int(max_diff_perc * len(string1))

        n_diff = 0
        for c1, c2 in zip(string1, string2):
            if c1 != c2:
                n_diff += 1
                if n_diff >= max_diff:
                    return False

        return True

    def __nearlyIdenticalGenomes(self, seqs, outputDir):
        identical = []
        numTaxa = 0

        nearlyIdenticalFile = os.path.join(outputDir, 'nearly_identical.tsv')
        if os.path.exists(nearlyIdenticalFile):
            for line in open(nearlyIdenticalFile):
                lineSplit = line.split('\t')
                s = set()
                for genomeId in lineSplit:
                    numTaxa += 1
                    s.add(genomeId.strip())
                identical.append(s)
        else:
            seqIds = list(seqs.keys())

            processed = set()
            for i in range(0, len(seqIds)):
                print('  %d of %d' % (i, len(seqIds)))
                seqIdI = seqIds[i]
                seqI = seqs[seqIdI]

                if seqIdI in processed:
                    continue

                processed.add(seqIdI)

                numTaxa += 1

                s = set()
                s.add(seqIdI)
                for j in range(i + 1, len(seqIds)):
                    seqIdJ = seqIds[j]
                    seqJ = seqs[seqIdJ]

                    if seqIdJ in processed:
                        continue

                    if self.__nearlyIdentical(seqI, seqJ):
                        s.add(seqIdJ)
                        processed.add(seqIdJ)

                identical.append(s)
                print('    set size: %d' % len(s))
                if len(s) > 1:
                    for genomeId in s:
                        genomeId = genomeId.replace('IMG_', '')
                        print(genomeId, self.metadata[genomeId]['taxonomy'])

            fout = open(nearlyIdenticalFile, 'w')
            for s in identical:
                fout.write('\t'.join(list(s)) + '\n')
            fout.close()

        print('  Number of taxa: %d' % numTaxa)
        print('  Number of dereplicated taxa: %d' % len(identical))

        return identical

    def run(self, outputDir):
        # make sure output directory exists
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)

        # remove similar taxa
        print('Filtering out highly similar taxa in order to reduce size of tree:')
        seqs = readFasta(self.derepConcatenatedAlignFile)

        nearlyIdentical = self.__nearlyIdenticalGenomes(seqs, outputDir)

        reducedSeqs = {}
        for s in nearlyIdentical:
            rndGenome = random.choice(tuple(s))
            reducedSeqs[rndGenome] = seqs[rndGenome]

        # write out reduced alignment
        reducedAlignmentFile = os.path.join(outputDir, "genome_tree.fasta")
        writeFasta(reducedSeqs, reducedAlignmentFile)

        # prune tree to retained taxa
        print('')
        print('Pruning tree:')
        tree = dendropy.Tree.get_from_path(self.tree, schema='newick', as_rooted=False, preserve_underscores=True)

        for seqId in reducedSeqs:
            node = tree.find_node_with_taxon_label(seqId)
            if not node:
                print('Missing taxa: %s' % seqId)

        tree.retain_taxa_with_labels(list(reducedSeqs.keys()))

        outputTree = os.path.join(outputDir, 'genome_tree.tre')
        tree.write_to_path(outputTree, schema='newick', suppress_rooting=True, unquoted_underscores=True)

        for t in tree.internal_nodes():
            t.label = None

        for t in tree.leaf_nodes():
            if t.taxon.label not in reducedSeqs:
                print('missing in sequence file: %s' % t.taxon.label)

        outputTreeWithoutLabels = os.path.join(outputDir, 'genome_tree.small.no_internal_labels.tre')
        tree.write_to_path(outputTreeWithoutLabels, schema='newick', suppress_rooting=True, unquoted_underscores=True)
        print('  Pruned tree written to: %s' % outputTree)

        # calculate model parameters for pruned tree
        print('')
        print('Determining model parameters for new tree.')
        outputTreeLog = os.path.join(outputDir, 'genome_tree.log')
        fastTreeOutput = os.path.join(outputDir, 'genome_tree.no_internal_labels.fasttree.tre')
        # os.system('FastTreeMP -nome -mllen -intree %s -log %s < %s > %s' % (outputTreeWithoutLabels, outputTreeLog, reducedAlignmentFile, fastTreeOutput))

        # calculate reference package for pruned tree
        print('')
        print('Creating reference package.')
        os.system('taxit create -l %s -P %s --aln-fasta %s --tree-stats %s --tree-file %s' % ('genome_tree_reduced', os.path.join(outputDir, 'genome_tree_reduced.refpkg'), reducedAlignmentFile, outputTreeLog, outputTree))

if __name__ == '__main__':
    print('RerootTree v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_tree_dir', help='genome tree directory for full tree')
    parser.add_argument('output_dir', help='output directory for reduced tree')

    args = parser.parse_args()

    createSamllTree = CreateSamllTree(args.genome_tree_dir)
    createSamllTree.run(args.output_dir)
