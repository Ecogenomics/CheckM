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

from lib.img import IMG

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
			genomeId, geneId = seqId.split('_')
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

	def run(self, geneTreeDir, extension, outputTree):
		# get genes in genomes
		print 'Reading all PFAM and TIGRFAM hits in trusted genomes.'
		img = IMG()
		genomeIds = img.genomeIds(img.genomeMetadata(), 'Trusted')
		genesInGenomes = self.__genesInGenomes(genomeIds)

		# read alignment files
		print 'Reading alignment files.'
		alignments = {}
		genomeIds = set()
		files = os.listdir(geneTreeDir)
		for f in files:
			if f.endswith(extension) and f != 'concatenated.faa':
				markerId = f[0:f.find('.')]
				print '  ' + f
				seqs = readFasta(os.path.join(geneTreeDir, f))
				seqs = self.__filterParalogs(seqs, markerId, genesInGenomes)

				genomeIds.update(set(seqs.keys()))
				alignments[f] = seqs

		# create concatenated alignment
		print 'Concatenating alignments.'
		concatenatedSeqs = {}
		for _, seqs in alignments.iteritems():
			alignLen = len(seqs[seqs.keys()[0]])
			for genomeId in genomeIds:
				if genomeId in seqs:
					# append alignment
					concatenatedSeqs[genomeId] = concatenatedSeqs.get(genomeId, '') + seqs[genomeId]
				else:
					# missing gene
					concatenatedSeqs[genomeId] = concatenatedSeqs.get(genomeId, '') + '-'*alignLen

		# save concatenated alignment
		concatenatedAlignFile = os.path.join(geneTreeDir, 'concatenated.faa')
		writeFasta(concatenatedSeqs, concatenatedAlignFile)

		# infer genome tree
		print 'Infering genome tree.'
		outputLog = outputTree[0:outputTree.rfind('.')] + '.log'
		cmd = 'FastTreeMP -wag -gamma -log ' + outputLog + ' ' + concatenatedAlignFile + ' > ' + outputTree
		os.system(cmd)

if __name__ == '__main__':
	print 'InferGenomeTree v' + __version__ + ': ' + __prog_desc__
	print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('gene_tree_dir', help='directory containing aligned genes to use for phylogenetic inference')
	parser.add_argument('-x', '--extension', help='extension of alignment files to process', default = '.faa')
	parser.add_argument('output_tree', help='output genome tree')

	args = parser.parse_args()

	inferGenomeTree = InferGenomeTree()
	inferGenomeTree.run(args.gene_tree_dir, args.extension, args.output_tree)
