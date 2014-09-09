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

"""
Find complete genomes in IMG.
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

class IdentifyCompleteGenomes(object):
    def __init__(self):
        pass

    def run(self, ubiquityThreshold, singleCopyThreshold, trustedCompleteness, trustedContamination, genomeCompleteness, genomeContamination):
        img = IMG()
        markerset = MarkerSet()

        metadata = img.genomeMetadata()

        trustedOut = open('./data/trusted_genomes.tsv', 'w')
        trustedOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tBiotic Relationship\tStatus\tCompleteness\tContamination\n')

        filteredOut = open('./data/filtered_genomes.tsv', 'w')
        filteredOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tBiotic Relationship\tStatus\tCompleteness\tContamination\n')

        allGenomeIds = set()
        allTrustedGenomeIds = set()
        for lineage in ['Archaea', 'Bacteria']:
            # get all genomes in lineage and build gene count table
            print '\nBuilding gene count table.'
            allLineageGenomeIds = img.genomeIdsByTaxonomy(lineage, metadata, 'All')
            countTable = img.countTable(allLineageGenomeIds)
            countTable = img.filterTable(allLineageGenomeIds, countTable, 0.9*ubiquityThreshold, 0.9*singleCopyThreshold)

            # get all genomes from specific lineage
            allGenomeIds = allGenomeIds.union(allLineageGenomeIds)

            print 'Lineage ' + lineage + ' contains ' + str(len(allLineageGenomeIds)) + ' genomes.'

            # tabulate genomes from each phylum
            allPhylumCounts = {}
            for genomeId in allLineageGenomeIds:
                taxon = metadata[genomeId]['taxonomy'][1]
                allPhylumCounts[taxon] = allPhylumCounts.get(taxon, 0) + 1

            # identify marker set for genomes
            markerGenes = markerset.markerGenes(allLineageGenomeIds, countTable, ubiquityThreshold*len(allLineageGenomeIds), singleCopyThreshold*len(allLineageGenomeIds))
            print '  Marker genes: ' + str(len(markerGenes))

            geneDistTable = img.geneDistTable(allLineageGenomeIds, markerGenes, spacingBetweenContigs=1e6)
            colocatedGenes = markerset.colocatedGenes(geneDistTable, metadata)
            colocatedSets = markerset.colocatedSets(colocatedGenes, markerGenes)
            print '  Marker set size: ' + str(len(colocatedSets))

            # identifying trusted genomes (highly complete, low contamination genomes)
            trustedGenomeIds = set()
            for genomeId in allLineageGenomeIds:
                completeness, contamination = markerset.genomeCheck(colocatedSets, genomeId, countTable)

                if completeness >= trustedCompleteness and contamination <= trustedContamination:
                    trustedGenomeIds.add(genomeId)
                    allTrustedGenomeIds.add(genomeId)

                    trustedOut.write(genomeId + '\t' + '; '.join(metadata[genomeId]['taxonomy']))
                    trustedOut.write('\t%.2f' % (float(metadata[genomeId]['genome size']) / 1e6))
                    trustedOut.write('\t' + str(metadata[genomeId]['scaffold count']))
                    trustedOut.write('\t' + metadata[genomeId]['biotic relationships'])
                    trustedOut.write('\t' + metadata[genomeId]['status'])
                    trustedOut.write('\t%.3f\t%.3f' % (completeness, contamination) + '\n')
                else:
                    filteredOut.write(genomeId + '\t' + '; '.join(metadata[genomeId]['taxonomy']))
                    filteredOut.write('\t%.2f' % (float(metadata[genomeId]['genome size']) / 1e6))
                    filteredOut.write('\t' + str(metadata[genomeId]['scaffold count']))
                    filteredOut.write('\t' + metadata[genomeId]['biotic relationships'])
                    filteredOut.write('\t' + metadata[genomeId]['status'])
                    filteredOut.write('\t%.3f\t%.3f' % (completeness, contamination) + '\n')

            print '  Trusted genomes: ' + str(len(trustedGenomeIds))

            # determine status of trusted genomes
            statusBreakdown = {}
            for genomeId in trustedGenomeIds:
                statusBreakdown[metadata[genomeId]['status']] = statusBreakdown.get(metadata[genomeId]['status'], 0) + 1

            print '  Trusted genome status breakdown: '
            for status, count in statusBreakdown.iteritems():
                print '    ' + status + ': ' + str(count)

            # determine status of retained genomes
            proposalNameBreakdown = {}
            for genomeId in trustedGenomeIds:
                proposalNameBreakdown[metadata[genomeId]['proposal name']] = proposalNameBreakdown.get(metadata[genomeId]['proposal name'], 0) + 1

            print '  Retained genome proposal name breakdown: '
            for pn, count in proposalNameBreakdown.iteritems():
                if 'KMG' in pn or 'GEBA' in pn or 'HMP' in pn:
                    print '    ' + pn + ': ' + str(count)

            print '  Filtered genomes by phylum:'
            trustedPhylumCounts = {}
            for genomeId in trustedGenomeIds:
                taxon = metadata[genomeId]['taxonomy'][1]
                trustedPhylumCounts[taxon] = trustedPhylumCounts.get(taxon, 0) + 1

            for phylum, count in allPhylumCounts.iteritems():
                print phylum + ': %d of %d' % (trustedPhylumCounts.get(phylum, 0), count)

        trustedOut.close()
        filteredOut.close()

        # write out lineage statistics for genome distribution
        allStats = {}
        trustedStats = {}

        for r in xrange(0, 6): # Domain to Genus
            for genomeId, data in metadata.iteritems():
                taxaStr = '; '.join(data['taxonomy'][0:r+1])
                allStats[taxaStr] = allStats.get(taxaStr, 0) + 1
                if genomeId in allTrustedGenomeIds:
                    trustedStats[taxaStr] = trustedStats.get(taxaStr, 0) + 1

        sortedLineages = img.lineagesSorted()

        fout = open('./data/lineage_stats.tsv', 'w')
        fout.write('Lineage\tGenomes with metadata\tTrusted genomes\n')
        for lineage in sortedLineages:
            fout.write(lineage + '\t' + str(allStats.get(lineage, 0))+ '\t' + str(trustedStats.get(lineage, 0))+ '\n')
        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find complete genomes in IMG.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--ubiquity', help='Ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='Single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-c', '--trusted_completeness', help='Completeness threshold for defining trusted genomes', type=float, default = 0.97)
    parser.add_argument('-d', '--trusted_contamination', help='Contamination threshold for defining trusted genomes', type=float, default = 0.03)
    parser.add_argument('-x', '--genome_completeness', help='Completeness threshold for defining final genome set', type=float, default = 0.95)
    parser.add_argument('-y', '--genome_contamination', help='Contamination threshold for defining final genome set', type=float, default = 0.05)

    args = parser.parse_args()

    identifyCompleteGenomes = IdentifyCompleteGenomes()
    identifyCompleteGenomes.run(args.ubiquity, args.single_copy, args.trusted_completeness, args.trusted_contamination, args.genome_completeness, args.genome_contamination)
