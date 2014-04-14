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

__prog_desc__ = 'infer genome tree from concatenation of genes'

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
from collections import defaultdict

from checkm.seqUtils import readFasta, writeFasta

from checkm.lib.img import IMG

from checkm.lib.taxonomyUtils import taxonomyWithRanks

class InferGenomeTree(object):
    def __init__(self):
        pass

    def __genesInGenomes(self, genomeIds):
        """ Get PFAM and TIGRFAM genes in genomes along with e-value. """
        genesInGenomes = {}
        for genomeId in genomeIds:
            markerIdToGeneIds = defaultdict(dict)
            with open(os.path.join(IMG.genomeDir, genomeId, genomeId + IMG.pfamExtension)) as f:
                f.readline()
                for line in f:
                    lineSplit = line.split('\t')
                    pfamId = lineSplit[8]
                    geneId = lineSplit[0]
                    evalue = float(lineSplit[6])
                    markerIdToGeneIds[pfamId][geneId] = min(evalue, markerIdToGeneIds[pfamId].get(geneId, 1000))

            with open(os.path.join(IMG.genomeDir, genomeId, genomeId + IMG.tigrExtension)) as f:
                f.readline()
                for line in f:
                    lineSplit = line.split('\t')
                    tigrId = lineSplit[6]
                    geneId = lineSplit[0]
                    evalue = float(lineSplit[4])
                    markerIdToGeneIds[tigrId][geneId] = min(evalue, markerIdToGeneIds[tigrId].get(geneId, 1000))

            genesInGenomes[genomeId] = markerIdToGeneIds

        return genesInGenomes

    def __filterParalogs(self, seqs, markerId, genesInGenomes):
        """ Identify paralogs and retain hit with lowest e-value."""

        # determine genome id and gene id of each sequences
        genomeIdGeneId = []
        for seqId in seqs.keys():
            genomeId, geneId = seqId.split('|')
            genomeIdGeneId.append([genomeId, geneId])

        # retain only best hit for paralogous genes
        filteredSeqs = {}
        for i, seqIdI in enumerate(seqs.keys()):
            genomeIdI, geneIdI = genomeIdGeneId[i]

            if genomeIdI in filteredSeqs:
                # this genome has already been processed
                continue

            geneIdToEvalue = genesInGenomes[genomeIdI][markerId]
            topHitId = seqIdI
            topHitEvalue = geneIdToEvalue[geneIdI]

            for j, seqIdJ in enumerate(seqs.keys()[i+1:]):
                genomeIdJ, geneIdJ = genomeIdGeneId[j]
                if genomeIdI == genomeIdJ:
                    if geneIdToEvalue[geneIdJ] < topHitEvalue:
                        topHitId = seqIdJ
                        topHitEvalue = geneIdToEvalue[geneIdJ]

            filteredSeqs[genomeIdI] = seqs[topHitId]

        return filteredSeqs
    
    def __taxonomy(self, img, genomeIds, outputTaxonomy):
        metadata = img.genomeMetadata()
        
        fout = open(outputTaxonomy, 'w')
        for genomeId in genomeIds:
            fout.write('IMG_' + genomeId + '\t')
            taxonomy = taxonomyWithRanks(metadata[genomeId]['taxonomy'])
            taxonomy = [x.replace('unclassified', '') for x in taxonomy]
            fout.write('; '.join(taxonomy) + '\n')
        
        fout.close()

    def run(self, geneTreeDir, alignmentDir, extension, outputAlignFile, outputTree, outputTaxonomy, bSupportValues = False):
        # read gene trees
        print 'Reading gene trees.'
        geneIds = set()
        files = os.listdir(geneTreeDir)
        for f in files:
            if f.endswith('.tre'):
                geneId = f[0:f.find('.')]
                geneIds.add(geneId)
                
        # write out genome tree taxonomy
        print 'Reading trusted genomes.'
        img = IMG()
        genomeIds = img.genomeMetadata().keys()
        self.__taxonomy(img, genomeIds, outputTaxonomy)
        
        print '  There are %d trusted genomes.' % (len(genomeIds))
    
        # get genes in genomes
        print 'Reading all PFAM and TIGRFAM hits in trusted genomes.'
        genesInGenomes = self.__genesInGenomes(genomeIds)

        # read alignment files
        print 'Reading alignment files.'
        alignments = {}
        genomeIds = set()
        files = os.listdir(alignmentDir)
        for f in files:
            geneId = f[0:f.find('.')]
            if f.endswith(extension) and geneId in geneIds:
                seqs = readFasta(os.path.join(alignmentDir, f))
                
                imgGeneId = geneId
                if imgGeneId.startswith('PF'):
                    imgGeneId = imgGeneId.replace('PF', 'pfam')
                seqs = self.__filterParalogs(seqs, imgGeneId, genesInGenomes)

                genomeIds.update(set(seqs.keys()))
                alignments[geneId] = seqs

        # create concatenated alignment
        print 'Concatenating alignments:'
        concatenatedSeqs = {}
        totalAlignLen = 0
        for geneId in sorted(alignments.keys()):
            seqs = alignments[geneId]
            alignLen = len(seqs[seqs.keys()[0]])
            print '  ' + str(geneId) + ',' + str(alignLen)
            totalAlignLen += alignLen
            for genomeId in genomeIds:
                if genomeId in seqs:
                    # append alignment
                    concatenatedSeqs['IMG_' + genomeId] = concatenatedSeqs.get('IMG_' + genomeId, '') + seqs[genomeId]
                else:
                    # missing gene
                    concatenatedSeqs['IMG_' + genomeId] = concatenatedSeqs.get('IMG_' + genomeId, '') + '-'*alignLen
                    
        print '  Total alignment length: ' + str(totalAlignLen)
        
        # save concatenated alignment
        writeFasta(concatenatedSeqs, outputAlignFile)

        # infer genome tree
        print 'Inferring genome tree.'
        outputLog = outputTree[0:outputTree.rfind('.')] + '.log'
        
        supportStr = ' '
        if not bSupportValues:
            supportStr = ' -nosupport '
        
        cmd = 'FastTreeMP' + supportStr + '-wag -gamma -log ' + outputLog + ' ' + outputAlignFile + ' > ' + outputTree
        os.system(cmd)

if __name__ == '__main__':
    print 'InferGenomeTree v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gene_tree_dir', help='directory containing gene trees to use for phylogenetic inference')
    parser.add_argument('alignment_dir', help='directory containing multiple sequence alignments for gene trees')
    parser.add_argument('output_align', help='output concatenated alignment')
    parser.add_argument('output_tree', help='output genome tree')
    parser.add_argument('output_taxonomy', help='output genome tree taxonomy')
    parser.add_argument('-x', '--extension', help='extension of alignment files to process', default = '.aln.masked.faa')

    args = parser.parse_args()

    inferGenomeTree = InferGenomeTree()
    inferGenomeTree.run(args.gene_tree_dir, args.alignment_dir, args.extension, args.output_align, args.output_tree, args.output_taxonomy)
