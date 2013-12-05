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

from lib.img import IMG
from lib.markerSet import MarkerSet

class QcGenomes(object):
    def __init__(self):
        pass
    
    def run(self, maxContigs, minN50, ubiquityThreshold, singleCopyThreshold, trustedCompleteness, trustedContamination):   
        img = IMG()
        markerset = MarkerSet()
        
        metadata = img.genomeMetadata()
          
        trustedOut = open('./data/trusted_genomes.tsv', 'w')
        trustedOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tBiotic Relationship\tStatus\tCompleteness\tContamination\n')
        
        filteredOut = open('./data/filtered_genomes.tsv', 'w')
        filteredOut.write('Genome Id\tLineage\tGenome size (Mbps)\tScaffold count\tBiotic Relationship\tStatus\tCompleteness\tContamination\n')
        
        allTrustedGenomeIds = set()
        for lineage in ['Archaea']: #, 'Bacteria']:  
            # get all genomes in lineage
            print 'Identifying all ' + lineage + ' genomes.'
            allLineageGenomeIds = img.genomeIdsByTaxonomy(lineage, metadata, 'All')
            print '  Number of genomes: ' + str(len(allLineageGenomeIds))
            
            # tabulate genomes from each phylum
            allPhylumCounts = {}
            for genomeId in allLineageGenomeIds:
                taxon = metadata[genomeId]['taxonomy'][1]
                allPhylumCounts[taxon] = allPhylumCounts.get(taxon, 0) + 1
            
            # filter out genomes missing sequence data 
            print '\nFilter out genomes missing data.'
            missing = img.genomesWithMissingData(allLineageGenomeIds)   
            filteredStatus = {}
            for genomeId in missing:
                filteredStatus[metadata[genomeId]['status']] = filteredStatus.get(metadata[genomeId]['status'], 0) + 1 
             
            print '  Filtered genomes: %d (%.2f%%)' % (len(missing), len(missing)*100.0 / len(allLineageGenomeIds))
            print '  ' + str(filteredStatus)
            retainedGenomeIds = allLineageGenomeIds - missing
            
            # filter out low-quality genomes
            print '\nFiltering out low quality genomes.'
            filteredGenomeIds = set()
            filteredStatus = {}
            for genomeId in retainedGenomeIds:
                if metadata[genomeId]['scaffold count'] > maxContigs or metadata[genomeId]['N50'] < minN50:
                    filteredStatus[metadata[genomeId]['status']] = filteredStatus.get(metadata[genomeId]['status'], 0) + 1 
                    filteredGenomeIds.add(genomeId)
                    if metadata[genomeId]['status'] == 'Finished':
                        print genomeId
                        print 'Scaffold count: ', metadata[genomeId]['scaffold count']
                        print 'N50: ', metadata[genomeId]['N50']
            
            print '  Filtered genomes: %d (%.2f%%)' % (len(filteredGenomeIds), len(filteredGenomeIds)*100.0 / len(allLineageGenomeIds))
            print '  ' + str(filteredStatus)
            retainedGenomeIds = retainedGenomeIds - filteredGenomeIds
                   
            # build gene count table
            print '\nDetermining initial marker gene sets for genome filtering.'
            countTable = img.countTable(retainedGenomeIds)
            countTable = img.filterTable(retainedGenomeIds, countTable, 0.9*ubiquityThreshold, 0.9*singleCopyThreshold)
                          
            # identify marker genes for genomes
            markerGenes = markerset.markerGenes(retainedGenomeIds, countTable, ubiquityThreshold*len(retainedGenomeIds), singleCopyThreshold*len(retainedGenomeIds))
            print '  Marker genes: ' + str(len(markerGenes))
                        
            tigrToRemove = img.identifyRedundantTIGRFAMs(markerGenes)
            markerGenes = markerGenes - tigrToRemove
            print '  Marker genes after filtering redundant TIGRFAMs: ' + str(len(markerGenes))
            
            fout = open('./data/trusted_marker_genes_' + lineage + '.txt', 'w')
            fout.write(str(markerGenes))
            fout.close()
            
            # identify marker sets
            geneDistTable = img.geneDistTable(retainedGenomeIds, markerGenes)
            colocatedGenes = markerset.colocatedGenes(geneDistTable)
            colocatedSets = markerset.colocatedSets(colocatedGenes, markerGenes)
            print '  Marker set size: ' + str(len(colocatedSets))
            
            fout = open('./data/trusted_marker_sets_' + lineage + '.txt', 'w')
            fout.write(str(colocatedSets))
            fout.close()
            
            # identifying trusted genomes (highly complete, low contamination genomes)
            print '\nIdentifying highly complete, low contamination genomes.'
            trustedGenomeIds = set()
            filteredGenomes = set()
            retainedStatus = {}
            filteredStatus = {}
            for genomeId in retainedGenomeIds:
                completeness, contamination = markerset.genomeCheck(colocatedSets, genomeId, countTable)
                 
                genomeStr = genomeId + '\t' + '; '.join(metadata[genomeId]['taxonomy'])
                genomeStr += '\t%.2f' % (float(metadata[genomeId]['genome size']) / 1e6)
                genomeStr +='\t' + str(metadata[genomeId]['scaffold count'])
                genomeStr +='\t' + metadata[genomeId]['biotic relationships']
                genomeStr +='\t' + metadata[genomeId]['status']
                genomeStr +='\t%.3f\t%.3f' % (completeness, contamination) + '\n'
                         
                if completeness >= trustedCompleteness and contamination <= trustedContamination:
                    trustedGenomeIds.add(genomeId)
                    allTrustedGenomeIds.add(genomeId)
                    retainedStatus[metadata[genomeId]['status']] = retainedStatus.get(metadata[genomeId]['status'], 0) + 1
                    
                    trustedOut.write(genomeStr)
                else:
                    filteredGenomes.add(genomeId)
                    filteredStatus[metadata[genomeId]['status']] = filteredStatus.get(metadata[genomeId]['status'], 0) + 1
                    
                    filteredOut.write(genomeStr)
                
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
    
    parser.add_argument('-m', '--max_contigs', help='maximum number of contigs/scaffolds to retain genome', type=int, default = 100)
    parser.add_argument('-n', '--min_N50', help='minimum N50 for retaining genome', type=int, default = 100000)
    parser.add_argument('-u', '--ubiquity', help='ubiquity threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-s', '--single_copy', help='single-copy threshold for defining marker set', type=float, default = 0.97)
    parser.add_argument('-c', '--trusted_completeness', help='completeness threshold for defining trusted genomes', type=float, default = 0.97)
    parser.add_argument('-d', '--trusted_contamination', help='contamination threshold for defining trusted genomes', type=float, default = 0.03)
    
    args = parser.parse_args()
    
    qcGenomes = QcGenomes()
    qcGenomes.run(args.max_contigs, args.min_N50, args.ubiquity, args.single_copy, args.trusted_completeness, args.trusted_contamination)
