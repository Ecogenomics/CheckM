#!/usr/bin/env python3

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

"""
Identify genomes that differ substantially from the infer marker set.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse

from lib.img import IMG
from lib.markerSet import MarkerSet

class IdentifyDegenerateGenomes(object):
    def __init__(self):
        pass

    def run(self, ubiquityThreshold, singleCopyThreshold, minGenomes, mostSpecificRank, minMarkers, completenessThreshold, contaminationThreshold):
        print(('Ubiquity threshold: ' + str(ubiquityThreshold)))
        print(('Single-copy threshold: ' + str(singleCopyThreshold)))
        print(('Min. genomes: ' + str(minGenomes)))
        print(('Most specific taxonomic rank: ' + str(mostSpecificRank)))
        print(('Min markers: ' + str(minMarkers)))
        print(('Completeness threshold: ' + str(completenessThreshold)))
        print(('Contamination threshold: ' + str(contaminationThreshold)))

        img = IMG()
        markerset = MarkerSet()

        lineages = img.lineagesByCriteria(minGenomes, mostSpecificRank)

        degenerateGenomes = {}
        for lineage in lineages:
            genomeIds = img.genomeIdsByTaxonomy(lineage, 'Final')

            print('')
            print(('Lineage ' + lineage + ' contains ' + str(len(genomeIds)) + ' genomes.'))

            # get table of PFAMs and do some initial filtering to remove PFAMs that are
            # clearly not going to pass the ubiquity and single-copy thresholds
            countTable = img.countTable(genomeIds)
            countTable = img.filterTable(genomeIds, countTable, ubiquityThreshold*0.9, singleCopyThreshold*0.9)

            markerGenes = markerset.markerGenes(genomeIds, countTable, ubiquityThreshold*len(genomeIds), singleCopyThreshold*len(genomeIds))
            if len(markerGenes) < minMarkers:
                continue

            geneDistTable = img.geneDistTable(genomeIds, markerGenes, spacingBetweenContigs=1e6)
            colocatedGenes = markerset.colocatedGenes(geneDistTable)
            colocatedSets = markerset.colocatedSets(colocatedGenes, markerGenes)

            for genomeId in genomeIds:
                completeness, contamination = markerset.genomeCheck(colocatedSets, genomeId, countTable)

                if completeness < completenessThreshold or contamination > contaminationThreshold:
                    degenerateGenomes[genomeId] = degenerateGenomes.get(genomeId, []) + [[lineage.split(';')[-1].strip(), len(genomeIds), len(colocatedSets), completeness, contamination]]

        # write out degenerate genomes
        metadata = img.genomeMetadata('Final')

        fout = open('./data/degenerate_genomes.tsv', 'w')
        fout.write('Genome Id\tTaxonomy\tGenome Size (Gbps)\tScaffolds\tBiotic Relationships\tStatus\tLineage\t# genomes\tMarker set size\tCompleteness\tContamination\n')
        for genomeId, data in list(degenerateGenomes.items()):
            fout.write(genomeId + '\t' + '; '.join(metadata[genomeId]['taxonomy']) + '\t%.2f' % (float(metadata[genomeId]['genome size']) / 1e6) + '\t' + str(metadata[genomeId]['scaffold count']))
            fout.write('\t' + metadata[genomeId]['biotic relationships'] + '\t' + metadata[genomeId]['status'])

            for d in data:
                fout.write('\t' + d[0] + '\t' + str(d[1]) + '\t' + str(d[2]) + '\t%.3f\t%.3f' % (d[3], d[4]))
            fout.write('\n')

        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify genomes that differ substantially from the infer marker set.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-g', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 40)
    parser.add_argument('-r', '--most_specific_rank', help='Most specific rank to include in analysis', type=int, default= 4)
    parser.add_argument('-m', '--min_markers', help='Minimum markers required to include in analysis', type=int, default = 20)
    parser.add_argument('-x', '--genome_completeness', help='Completeness threshold for defining final genome set', type=float, default = 0.95)
    parser.add_argument('-y', '--genome_contamination', help='Contamination threshold for defining final genome set', type=float, default = 0.05)

    args = parser.parse_args()

    identifyDegenerateGenomes = IdentifyDegenerateGenomes()
    identifyDegenerateGenomes.run(args.ubiquity, args.single_copy, args.min_genomes, args.most_specific_rank, args.min_markers, args.genome_completeness, args.genome_contamination)
