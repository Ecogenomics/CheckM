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

__prog_desc__ = 'identify genomes suitable for calculating marker genes'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse
from collections import defaultdict

from checkm.lib.img import IMG
from lib.markerSetBuilder import MarkerSetBuilder

class QcGenomes(object):
    def __init__(self):
        pass

    def __genomeString(self, genomeId, metadata, completeness, contamination, missingMarkers, duplicateMarkers):
        genomeStr = genomeId + '\t' + '; '.join(metadata[genomeId]['taxonomy'])
        genomeStr += '\t%.2f' % (float(metadata[genomeId]['genome size']) / 1e6)
        genomeStr +='\t' + str(metadata[genomeId]['scaffold count'])
        genomeStr +='\t' + str(metadata[genomeId]['gene count'])
        genomeStr +='\t' + str(metadata[genomeId]['coding base count'])
        genomeStr +='\t' + str(metadata[genomeId]['N50'])
        genomeStr +='\t' + metadata[genomeId]['biotic relationships']
        genomeStr +='\t' + metadata[genomeId]['status']
        genomeStr +='\t%.3f\t%.3f' % (completeness, contamination)
        genomeStr += '\t' + ', '.join(list(missingMarkers))
        genomeStr += '\t' + ', '.join(list(duplicateMarkers)) + '\n'

        return genomeStr

    def run(self, inputMetadataFile, outputMetadataFile, ubiquityThreshold, singleCopyThreshold, trustedCompleteness, trustedContamination):
        img = IMG()
        markerSetBuilder = MarkerSetBuilder()

        allOut = open('./data/genomes_all.tsv', 'w')
        allOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tGene count\tCoding base count\tN50\tBiotic Relationship\tStatus\tCompleteness\tContamination\tMissing markers\tDuplicate markers\n')

        trustedOut = open('./data/genomes_trusted.tsv', 'w')
        trustedOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tGene count\tCoding base count\tN50\tBiotic Relationship\tStatus\tCompleteness\tContamination\tMissing markers\tDuplicate markers\n')

        filteredOut = open('./data/genomes_filtered.tsv', 'w')
        filteredOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tGene count\tCoding base count\tN50\tBiotic Relationship\tStatus\tCompleteness\tContamination\tMissing markers\tDuplicate markers\n')

        metadataOut = open(outputMetadataFile, 'w')
        
        # read input metadata file
        metadata = img.genomeMetadataFromFile(inputMetadataFile)
        
        finishedGenomes = defaultdict(set)
        allGenomes = defaultdict(set)
        
        metadataLine = {}
        
        bHeader = True
        for line in open(inputMetadataFile):
            if bHeader:
                metadataOut.write(line)
                bHeader = False
                continue
            
            lineSplit = line.split('\t')
            genomeId = lineSplit[0]
            domain = lineSplit[1]
            status = lineSplit[2]
            
            if status == 'Finished':
                finishedGenomes[domain].add(genomeId)
            
            allGenomes[domain].add(genomeId)
            metadataLine[genomeId] = line

        allTrustedGenomeIds = set()
        for lineage, allLineageGenomeIds in allGenomes.iteritems():
            print '[' + lineage + ']'
            print '  Number of genomes: %d' % len(allLineageGenomeIds)

            # tabulate genomes from each phylum
            allPhylumCounts = {}
            for genomeId in allLineageGenomeIds:
                taxon = metadata[genomeId]['taxonomy'][1]
                allPhylumCounts[taxon] = allPhylumCounts.get(taxon, 0) + 1

            # identify marker genes for finished genomes
            print '\nDetermining initial marker gene sets for genome filtering.'
            markerSet = markerSetBuilder.buildMarkerSet(finishedGenomes[lineage], ubiquityThreshold, singleCopyThreshold)

            print '  Marker set consists of %s marker genes organized into %d sets.' % (markerSet.numMarkers(), markerSet.numSets())
            fout = open('./data/trusted_marker_sets_' + lineage + '.txt', 'w')
            fout.write(str(markerSet.markerSet))
            fout.close()

            # identifying trusted genomes (highly complete, low contamination genomes)
            print '\nIdentifying highly complete, low contamination genomes.'
            trustedGenomeIds = set()
            filteredGenomes = set()
            retainedStatus = {}
            filteredStatus = {}
            geneCountTable = img.geneCountTable(allLineageGenomeIds)
            for genomeId in allLineageGenomeIds:
                completeness, contamination, missingMarkers, duplicateMarkers = markerSetBuilder.genomeCheck(markerSet.markerSet, genomeId, geneCountTable)
                
                genomeStr = self.__genomeString(genomeId, metadata, completeness, contamination, missingMarkers, duplicateMarkers)

                if completeness >= trustedCompleteness and contamination <= trustedContamination:
                    trustedGenomeIds.add(genomeId)
                    allTrustedGenomeIds.add(genomeId)
                    retainedStatus[metadata[genomeId]['status']] = retainedStatus.get(metadata[genomeId]['status'], 0) + 1

                    trustedOut.write(genomeStr)
                    allOut.write(genomeStr)
                    
                    metadataOut.write(metadataLine[genomeId])
                else:
                    filteredGenomes.add(genomeId)
                    filteredStatus[metadata[genomeId]['status']] = filteredStatus.get(metadata[genomeId]['status'], 0) + 1

                    filteredOut.write(genomeStr)
                    allOut.write(genomeStr)

            print '  Filtered genomes: %d (%.2f%%)' % (len(filteredGenomes), len(filteredGenomes)*100.0 / len(allLineageGenomeIds))
            print '  ' + str(filteredStatus)
            print '  \nTrusted genomes: %d (%.2f%%)' % (len(trustedGenomeIds), len(trustedGenomeIds)*100.0 / len(allLineageGenomeIds))
            print '  ' + str(retainedStatus)

            # determine status of retained genomes
            print '\nTrusted genomes by phylum:'
            trustedPhylumCounts = {}
            for genomeId in trustedGenomeIds:
                taxon = metadata[genomeId]['taxonomy'][1]
                trustedPhylumCounts[taxon] = trustedPhylumCounts.get(taxon, 0) + 1

            for phylum, count in allPhylumCounts.iteritems():
                print '  ' + phylum + ': %d of %d' % (trustedPhylumCounts.get(phylum, 0), count)
            print ''

        allOut.close()
        trustedOut.close()
        filteredOut.close()
        metadataOut.close()

        # write out lineage statistics for genome distribution
        allStats = {}
        trustedStats = {}

        for r in xrange(0, 6): # Domain to Genus
            for genomeId, data in metadata.iteritems():
                taxaStr = ';'.join(data['taxonomy'][0:r+1])
                allStats[taxaStr] = allStats.get(taxaStr, 0) + 1
                if genomeId in allTrustedGenomeIds:
                    trustedStats[taxaStr] = trustedStats.get(taxaStr, 0) + 1

        sortedLineages = img.lineagesSorted(metadata)

        fout = open('./data/lineage_stats.tsv', 'w')
        fout.write('Lineage\tGenomes with metadata\tTrusted genomes\n')
        for lineage in sortedLineages:
            fout.write(lineage + '\t' + str(allStats.get(lineage, 0))+ '\t' + str(trustedStats.get(lineage, 0))+ '\n')
        fout.close()

if __name__ == '__main__':
    print 'QcGenomes v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', help='input IMG metadata file')
    parser.add_argument('output_file', help='new IMG metadata file')

    parser.add_argument('-u', '--ubiquity', help='ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-c', '--trusted_completeness', help='completeness threshold for defining trusted genomes', type=float, default = 0.97)
    parser.add_argument('-d', '--trusted_contamination', help='contamination threshold for defining trusted genomes', type=float, default = 0.03)

    args = parser.parse_args()

    qcGenomes = QcGenomes()
    qcGenomes.run(args.input_file, args.output_file, args.ubiquity, args.single_copy, args.trusted_completeness, args.trusted_contamination)
