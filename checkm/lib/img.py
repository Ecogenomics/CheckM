###############################################################################
#
# pfam.py - utilities for interfacing with IMG data.
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
import random

class IMG(object):
    metadataFile = '/srv/whitlam/bio/db/img/02122013/metadata/img_metadata.latest.tsv'
    genomeDir = '/srv/whitlam/bio/db/img/02122013/genomes/'
    pfamExtension = '.pfam.tab.txt'
    tigrExtension = '.tigrfam.tab.txt'
    pfamHMMs = '/srv/whitlam/bio/db/pfam/27/Pfam-A.hmm'
    tigrHMMs = '/srv/whitlam/bio/db/tigrfam/13.0/Tigr.hmm'

    def __init__(self):
        self.numScaffoldThreshold = 300
        self.N50Threshold = 10000
        pass

    def trustedGenomes(self):
        genomeIds = set()
        with open('./data/genomes_trusted.tsv') as f:
            f.readline()
            for line in f:
                genomeIds.add(line.split('\t')[0])

        return genomeIds

    def genomeIds(self, metadata, genomeFilter = 'All'):
        genomeIds = set()

        if genomeFilter == 'All':
            for genomeId, data in metadata.iteritems():
                if data['Donovan Note'] == '':
                    genomeIds.add(genomeId)
        elif genomeFilter == 'Trusted':
                genomeIds = self.trustedGenomes()
        else:
            print '[Error] Unrecognized filter flag: ' + genomeFilter
            sys.exit()

        return genomeIds

    def geneIdToScaffoldId(self, genomeId):
        d = {}

        for line in open(self.genomeDir + genomeId + '/' + genomeId + '.gff'):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')

            scaffoldId = lineSplit[0]
            geneId = lineSplit[8]
            geneId = geneId[geneId.find('=')+1:geneId.find(';')]

            d[geneId] = scaffoldId

        return d

    def pfamIdToGeneId(self, genomeId):
        return self.clusterIdToGeneId(genomeId, self.pfamExtension)

    def tigrIdToGeneId(self, genomeId):
        return self.clusterIdToGeneId(genomeId, self.tigrExtension)

    def clusterIdToGeneId(self, genomeId, extension):
        d = {}

        bHeader = True
        for line in open(self.genomeDir + genomeId + '/' + genomeId + extension):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            geneId = lineSplit[0]
            clusterId = lineSplit[8]

            d[clusterId] = d.get(clusterId, []) + [geneId]

        return d

    def genomeIdToTaxonomy(self):
        metadata = self.genomeMetadata()

        d= {}
        for genomeId in metadata:
            d[genomeId] = '; '.join(metadata[genomeId]['taxonomy'])

        return d

    def genomeMetadata(self):
        metadata = {}

        bHeader = True
        for line in open(IMG.metadataFile):
            lineSplit = line.split('\t')
            lineSplit = [x.strip() for x in lineSplit]

            if bHeader:
                statusIndex = lineSplit.index('Status')
                proposalNameIndex = lineSplit.index('Proposal Name')
                scaffoldCountIndex = lineSplit.index('Scaffold Count')
                genomeSizeIndex = lineSplit.index('Genome Size')
                geneCountIndex = lineSplit.index('Gene Count')
                codingBaseCountIndex = lineSplit.index('Coding Base Count')
                bioticRelationshipsIndex = lineSplit.index('Biotic Relationships')
                n50Index = lineSplit.index('N50')
                donovanNoteIndex = lineSplit.index('Donovan Note')
                bHeader = False
                continue

            genomeId = lineSplit[0].strip()

            rDomain = lineSplit[1].strip()
            rPhylum = lineSplit[6].strip()
            rClass = lineSplit[7].strip()
            rOrder = lineSplit[8].strip()
            rFamily = lineSplit[9].strip()
            rGenus = lineSplit[10].strip()
            rSpecies = lineSplit[11].strip()

            metadata[genomeId] = {}
            metadata[genomeId]['status'] = lineSplit[statusIndex]
            metadata[genomeId]['proposal name'] = lineSplit[proposalNameIndex]
            metadata[genomeId]['taxonomy'] = [rDomain, rPhylum, rClass, rOrder, rFamily, rGenus, rSpecies]
            metadata[genomeId]['scaffold count'] = int(lineSplit[scaffoldCountIndex])

            try:
                metadata[genomeId]['genome size'] = int(lineSplit[genomeSizeIndex])
            except:
                metadata[genomeId]['genome size'] = 'NA'

            try:
                metadata[genomeId]['gene count'] = int(lineSplit[geneCountIndex])
            except:
                metadata[genomeId]['gene count'] = 'NA'

            try:
                metadata[genomeId]['coding base count'] = int(lineSplit[codingBaseCountIndex])
            except:
                metadata[genomeId]['coding base count'] = 'NA'

            metadata[genomeId]['biotic relationships'] = lineSplit[bioticRelationshipsIndex]
            metadata[genomeId]['N50'] = int(lineSplit[n50Index])
            metadata[genomeId]['Donovan Note'] = lineSplit[donovanNoteIndex]

        return metadata

    def genomesWithMissingData(self, genomeIds):
        missingPFAM = self.missingPfamData(genomeIds)
        missingTIGR = self.missingTigrData(genomeIds)

        return missingPFAM.union(missingTIGR)

    def missingPfamData(self, genomeIds):
        missing = set()
        for genomeId in genomeIds:
            if not os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + self.pfamExtension):
                missing.add(genomeId)

                #if os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + '.genes.fna'):
                #  print '[Warning] ' + genomeId + ' contains ORF data, but not PFAM annotations.'

        return missing

    def missingTigrData(self, genomeIds):
        missing = set()
        for genomeId in genomeIds:
            if not os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + self.tigrExtension):
                missing.add(genomeId)

                #if os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + '.genes.fna'):
                #  print '[Warning] ' + genomeId + ' contains ORF data, but not TIGRFAM annotations.'

        return missing

    def genomeIdsByTaxonomy(self, taxonStr, metadata, genomeFilter = 'All'):
        genomeIds = self.genomeIds(metadata, genomeFilter)

        searchTaxa = taxonStr.split(';')
        genomeIdsOfInterest = set()
        for genomeId in genomeIds:
            bKeep = True
            for r in xrange(0, len(searchTaxa)):
                if taxonStr == 'universal':
                    bKeep = True
                elif taxonStr == 'prokaryotes' and (metadata[genomeId]['taxonomy'][0] == 'Bacteria' or metadata[genomeId]['taxonomy'][0] == 'Archaea'):
                    bKeep = True
                elif searchTaxa[r].strip() == metadata[genomeId]['taxonomy'][r]:
                    bKeep = True
                else:
                    bKeep = False
                    break

            if bKeep:
                genomeIdsOfInterest.add(genomeId)

        return genomeIdsOfInterest

    def lineageStats(self):
        metadata = self.genomeMetadata()

        stats = {}
        for r in xrange(0, 6): # Domain to Genus
            for _, data in metadata.iteritems():
                taxaStr = '; '.join(data['taxonomy'][0:r+1])
                stats[taxaStr] = stats.get(taxaStr, 0) + 1

        return stats

    def lineagesSorted(self, metadata, mostSpecificRank=5):
        lineages = []
        for r in xrange(0, mostSpecificRank+1):
            taxa = set()
            for _, data in metadata.iteritems():
                if 'unclassified' not in data['taxonomy'][0:r+1]:
                    taxa.add('; '.join(data['taxonomy'][0:r+1]))

            lineages += sorted(list(taxa))

        return lineages

    def lineagesByCriteria(self, minGenomes, mostSpecificRank):
        l = []

        stats = self.lineageStats()
        for lineage in self.lineagesSorted(mostSpecificRank):
            if stats[lineage] > minGenomes:
                l.append(lineage)

        return l

    def __readTable(self, table, genomeIds, extension, clusterIdIndex):
        for genomeId in genomeIds:
            count = {}
            bHeader = True
            for line in open(os.path.join(self.genomeDir, genomeId, genomeId + extension)):
                if bHeader:
                    bHeader = False
                    continue

                lineSplit = line.split('\t')
                clusterId = lineSplit[clusterIdIndex]
                count[clusterId] = count.get(clusterId, 0) + 1

            for clusterId, c in count.iteritems():
                if clusterId not in table:
                    table[clusterId] = {}
                table[clusterId][genomeId] = c

    def countTable(self, genomeIds):
        table = {}
        self.__readTable(table, genomeIds, self.pfamExtension, 8)
        self.__readTable(table, genomeIds, self.tigrExtension, 6)

        return table

    def filterTable(self, genomeIds, table, ubiquityThreshold = 0.9, singleCopyThreshold = 0.9):
        idsToFilter = []
        for pfamId, genomeCounts in table.iteritems():
            ubiquity = 0
            singleCopy = 0
            for genomeId in genomeIds:
                count = genomeCounts.get(genomeId, 0)

                if count > 0:
                    ubiquity += 1

                if count == 1:
                    singleCopy += 1

            if (float(ubiquity)/len(genomeIds) < ubiquityThreshold) or (float(singleCopy)/len(genomeIds) < singleCopyThreshold):
                idsToFilter.append(pfamId)

        for clusterId in idsToFilter:
            table.pop(clusterId)

        return table

    def clusterIdToGeneIds(self, markerGenes, filename, clusterIdIndex):
        bHeader = True
        geneIdsToClusterIds = {}
        for line in open(filename):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            geneId = lineSplit[0]
            clusterId = lineSplit[clusterIdIndex]
            if clusterId in markerGenes:
                geneIdsToClusterIds[clusterId] = geneIdsToClusterIds.get(clusterId, []) + [geneId]

        return geneIdsToClusterIds

    def geneDistTable(self, genomeIds, markerGenes):
        table = {}
        for genomeId in genomeIds:
            # read cluster table for gene identifiers
            filename = self.genomeDir + genomeId + '/' + genomeId
            pfamIdToGeneIds = self.clusterIdToGeneIds(markerGenes, filename + self.pfamExtension, 8)
            tigrIdToGeneIds = self.clusterIdToGeneIds(markerGenes, filename + self.tigrExtension, 6)

            # read gene locations
            genePosition = {}
            gffFile = self.genomeDir + genomeId + '/' + genomeId + '.gff'
            for line in open(gffFile):
                if line[0] == '#':
                    continue

                lineSplit = line.split('\t')

                geneId = lineSplit[8].split(';')[0]
                geneId = geneId[geneId.find('=')+1:]

                start = int(lineSplit[3])
                end = int(lineSplit[4])

                genePosition[geneId] = [start, end]

            # create gene mapping table
            clusterIdToGenomePositions = {}
            for pfamId, geneIds in pfamIdToGeneIds.iteritems():
                positions = []
                for geneId in geneIds:
                    positions += [genePosition[geneId]]
                clusterIdToGenomePositions[pfamId] = positions

            for tigrId, geneIds in tigrIdToGeneIds.iteritems():
                positions = []
                for geneId in geneIds:
                    positions += [genePosition[geneId]]
                clusterIdToGenomePositions[tigrId] = positions

            table[genomeId] = clusterIdToGenomePositions

        return table

    def sampleGenome(self, genomeLen, percentComplete, contigLen):
        contigsInGenome = genomeLen / contigLen

        contigsToSample = int(contigsInGenome*percentComplete + 0.5)

        sampledContigs = random.sample(xrange(contigsInGenome), contigsToSample)
        sampledContigs.sort()

        for i in xrange(contigsToSample):
            sampledContigs[i] *= contigLen

        return sampledContigs

    def containedMarkerGenes(self, markerGenes, clusterIdToGenomePositions, startPartialGenomeContigs, contigLen):
        contained = set()

        for markerGene in markerGenes:
            positions = clusterIdToGenomePositions.get(markerGene, [])

            bContained = False
            for p in positions:
                for s in startPartialGenomeContigs:
                    if (p[0] - s) >= 0 and (p[0] - s) < contigLen:
                        bContained = True
                        break

                if bContained:
                    break

            if bContained:
                contained.add(markerGene)

        return contained

    def identifyMitochondrialChloroplastGenes(self, genomeId):
        # identify mitochondrial or chloroplast sequences
        mitoChloroSeqs = set()
        for line in open(self.genomeDir + genomeId + '/' + genomeId + '.fna'):
            if line[0] == '>':
                if 'mitochondria' in line.lower() or 'mitochondrion' in line.lower() or 'chloroplast' in line.lower():
                    mitoChloroSeqs.add(line[1:].split()[0])

        # identify mitochondrial or chloroplast genes
        mitoChloroGenes = set()
        for line in open(self.genomeDir + genomeId + '/' + genomeId + '.gff'):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId in mitoChloroSeqs:
                desc = lineSplit[8]
                geneId = desc.split(';')[0].split('=')[1]
                mitoChloroGenes.add(geneId)

        return mitoChloroGenes

    def identifyRedundantTIGRFAMs(self, markerGenes):
        tigrIdToPfamId = {}
        for line in open('../data/tigrfam2pfam.tsv'):
            lineSplit = line.split('\t')
            pfamId = lineSplit[0]
            tigrId = lineSplit[1].rstrip()

            tigrIdToPfamId[tigrId] = tigrIdToPfamId.get(tigrId, []) + [pfamId]


        tigrToRemove = set()
        for markerGene in markerGenes:
            if markerGene in tigrIdToPfamId:
                for pfamId in tigrIdToPfamId[markerGene]:
                    if pfamId in markerGenes:
                        tigrToRemove.add(markerGene)

        return tigrToRemove
