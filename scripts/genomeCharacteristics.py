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
Calculate genome characteristics (coding density, strandedness, start codon frequency,
  stop codon frequency, gene length distributions) under translation table 11 for
  different groups of organisms (bacteria, archaea, bacteriophage, fungi, protists).
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os, sys
import string

from lib.img import IMG

from numpy import mean, std

class GenomeCharacteristics(object):
    def __init__(self):
        self.complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
        pass

    def revComp(self, dna):
        return dna.translate(self.complements)[::-1]

    def runProdigal(self, genomeId):
        print '  Running Prodigal.'
        #os.system('prodigal -c -f gff -g 11 -d ./genome_characteristics/' + genomeId + '.genes.fna -a ./genome_characteristics/' + genomeId + '.genes.faa -i ' + IMG.genomeDir + genomeId + '/' + genomeId + '.fna > ./genome_characteristics/' + genomeId + '.prodigal.gff\n')

    def readSeqs(self, filename):
        print '  Reading sequences.'
        seqs = {}
        seqId = None
        seq = ''
        for line in open(filename):
            if line[0] == '>':
                if seqId != None:
                    seqs[seqId] = seq

                seqId = line[1:].split()[0]
                seq = ''
            else:
                seq += line.rstrip()

        seqs[seqId] = seq

        return seqs

    def parseLineGFF(self, line, genes, posStrand, negStrand, startCodonFreq, stopCodonFreq, rbsMotifFreq, geneLens, seqs):
        lineSplit = line.split('\t')

        genes += 1

        seqId = lineSplit[0]
        seq = seqs[seqId]

        start = int(lineSplit[3])
        end = int(lineSplit[4])
        geneLen = end - start + 1
        geneLens.append(geneLen)

        strand = lineSplit[6]
        if strand == '+':
            posStrand += 1
            stopCodon = seq[end-3:end]
        elif strand == '-':
            negStrand += 1
            stopCodon = self.revComp(seq[start-1:start+2])
        else:
            print '[Error] Unknown strand type.'
            sys.exit()

        features = lineSplit[8].split(';')
        start_type = features[2].split('=')[1]
        startCodonFreq[start_type] = startCodonFreq.get(start_type, 0) + 1

        stopCodonFreq[stopCodon] = stopCodonFreq.get(stopCodon, 0) + 1

        rbs_type = features[3].split('=')[1]
        rbsMotifFreq[rbs_type] = rbsMotifFreq.get(rbs_type, 0) + 1

        return genes, posStrand, negStrand

    def parseGFF(self, genomeId):
        seqs = self.readSeqs(IMG.genomeDir + genomeId + '/' + genomeId + '.fna')

        genes = 0
        genesPlasmid = 0

        posStrand = 0
        posStrandPlasmid = 0

        negStrand = 0
        negStrandPlasmid = 0

        startCodonFreq = {}
        startCodonFreqPlasmid = {}

        stopCodonFreq = {}
        stopCodonFreqPlasmid = {}

        rbsMotifFreq = {}
        rbsMotifFreqPlasmid = {}

        geneLens = []
        geneLensPlasmid = []

        totalSeqLen = 0
        totalSeqLenPlasmid = 0

        bPlasmid = False
        lineNum = 0
        for line in open('./genome_characteristics/' + genomeId + '.prodigal.gff'):
            lineNum += 1

            if lineNum % 10000 == 0:
                print '  Reading line: ' + str(lineNum)

            if line[0] == '#':
                if 'seqlen=' in line:
                    bPlasmid = ('plasmid' in line.lower())
                    lineSplit = line.split(';')
                    seqLen = int(lineSplit[1].split('=')[1])

                    if bPlasmid:
                        totalSeqLenPlasmid += seqLen
                    else:
                        totalSeqLen += seqLen

                continue

            if bPlasmid:
                genesPlasmid, posStrandPlasmid, negStrandPlasmid = self.parseLineGFF(line, genesPlasmid, posStrandPlasmid, negStrandPlasmid, startCodonFreqPlasmid, stopCodonFreqPlasmid, rbsMotifFreqPlasmid, geneLensPlasmid, seqs)
            else:
                genes, posStrand, negStrand = self.parseLineGFF(line, genes, posStrand, negStrand, startCodonFreq, stopCodonFreq, rbsMotifFreq, geneLens, seqs)

        if totalSeqLen != 0 and genes == 0:
            print '[Error] No Genes!'
            sys.exit()

        codingDensity = float( sum(geneLens) ) / max(totalSeqLen, 1)
        codingDensityPlasmid = float( sum(geneLensPlasmid) ) / max(totalSeqLenPlasmid, 1)
        for codon, count in startCodonFreq.iteritems():
            startCodonFreq[codon] = float( count ) / max(genes, 1)
        for codon, count in stopCodonFreq.iteritems():
            stopCodonFreq[codon] = float( count ) / max(genes, 1)
        for motif, count in rbsMotifFreq.iteritems():
            rbsMotifFreq[motif] = float( count ) / max(genes, 1)

        strandedness = float( abs(posStrand - negStrand) ) / max(genes, 1)
        strandednessPlasmid = float( abs(posStrandPlasmid - negStrandPlasmid) ) / max(genesPlasmid, 1)
        for codon, count in startCodonFreqPlasmid.iteritems():
            startCodonFreqPlasmid[codon] = float( count ) / max(genesPlasmid, 1)
        for codon, count in stopCodonFreqPlasmid.iteritems():
            stopCodonFreqPlasmid[codon] = float( count ) / max(genesPlasmid, 1)
        for motif, count in rbsMotifFreqPlasmid.iteritems():
            rbsMotifFreqPlasmid[motif] = float( count ) / max(genes, 1)

        return [[totalSeqLen, codingDensity, strandedness, startCodonFreq, stopCodonFreq, rbsMotifFreq, geneLens], [totalSeqLenPlasmid, codingDensityPlasmid, strandednessPlasmid, startCodonFreqPlasmid, stopCodonFreqPlasmid, rbsMotifFreqPlasmid, geneLensPlasmid]];

    def getCharacteristics(self, genomeId, groupName, charByGroup):
        genomicChar, plasmidChar = self.parseGFF(genomeId)

        if genomicChar[0] == 0: # indicates file contained not sequence information or insufficent data to run prodigal
            return

        if groupName not in charByGroup:
            charByGroup[groupName] = {}

        charByGroup[groupName]['seq_len'] = charByGroup[groupName].get('seq_len', []) + [genomicChar[0]]
        charByGroup[groupName]['coding_density'] = charByGroup[groupName].get('coding_density', []) + [genomicChar[1]]
        charByGroup[groupName]['strandedness'] = charByGroup[groupName].get('strandedness', []) + [genomicChar[2]]
        charByGroup[groupName]['start_codon_freq'] = charByGroup[groupName].get('start_codon_freq', []) + [genomicChar[3]]
        charByGroup[groupName]['stop_codon_freq'] = charByGroup[groupName].get('stop_codon_freq', []) + [genomicChar[4]]
        charByGroup[groupName]['rbs_motif'] = charByGroup[groupName].get('rbs_motif', []) + [genomicChar[5]]
        charByGroup[groupName]['gene_lens'] = charByGroup[groupName].get('gene_lens', []) + [genomicChar[6]]

        if plasmidChar[0] != 0:
            plasmidGroupName = groupName + ' plasmids'
            if plasmidGroupName not in charByGroup:
                charByGroup[groupName + ' plasmids'] = {}

            charByGroup[plasmidGroupName]['seq_len'] = charByGroup[plasmidGroupName].get('seq_len', []) + [plasmidChar[0]]
            charByGroup[plasmidGroupName]['coding_density'] = charByGroup[plasmidGroupName].get('coding_density', []) + [plasmidChar[1]]
            charByGroup[plasmidGroupName]['strandedness'] = charByGroup[plasmidGroupName].get('strandedness', []) + [plasmidChar[2]]
            charByGroup[plasmidGroupName]['start_codon_freq'] = charByGroup[plasmidGroupName].get('start_codon_freq', []) + [plasmidChar[3]]
            charByGroup[plasmidGroupName]['stop_codon_freq'] = charByGroup[plasmidGroupName].get('stop_codon_freq', []) + [plasmidChar[4]]
            charByGroup[plasmidGroupName]['rbs_motif'] = charByGroup[plasmidGroupName].get('rbs_motif', []) + [plasmidChar[5]]
            charByGroup[plasmidGroupName]['gene_lens'] = charByGroup[plasmidGroupName].get('gene_lens', []) + [plasmidChar[6]]

    def writeResults(self, charByGroup):
        fout = open('./data/genome_characteristics.tsv', 'w')
        fout.write('Group name\t# genomes\tSeq length\t\tCoding density\t\tStrandedness\t\tGene length\t\tStart codon freq\t\t\t\t\t\tStop codon freq\t\tRBS motif freq\n')
        for groupName, data in charByGroup.iteritems():
            #sanity check
            if len(data['seq_len']) != len(data['coding_density']) or len(data['seq_len']) != len(data['strandedness']) or len(data['seq_len']) != len(data['gene_lens']):
                print '[Error] Not all genome characteristics have the same length: ' + groupName
                sys.exit()

            fout.write(groupName)
            fout.write('\t' + str(len(data['seq_len'])))
            fout.write('\t%.2f\t%.3f' % (mean(data['seq_len']), std(data['seq_len'])))
            fout.write('\t%.2f\t%.3f' % (mean(data['coding_density']), std(data['coding_density'])))
            fout.write('\t%.2f\t%.3f' % (mean(data['strandedness']), std(data['strandedness'])))

            meanGeneLen = []
            for geneLenList in data['gene_lens']:
                meanGeneLen.append(mean(geneLenList))

            fout.write('\t%.2f\t%.3f' % (mean(meanGeneLen), std(meanGeneLen)))

            meanStartCodonDic = {}
            for startCodonDic in data['start_codon_freq']:
                for codon, freq in startCodonDic.iteritems():
                    meanStartCodonDic[codon] = meanStartCodonDic.get(codon, []) + [freq]

            for codon, freqList in meanStartCodonDic.iteritems():
                fout.write('\t' + codon + ': %.2f\t%.3f' % (mean(freqList), std(freqList)))

            meanStopCodonDic = {}
            for stopCodonDic in data['stop_codon_freq']:
                for codon, freq in stopCodonDic.iteritems():
                    meanStopCodonDic[codon] = meanStopCodonDic.get(codon, []) + [freq]

            for codon, freqList in meanStopCodonDic.iteritems():
                fout.write('\t' + codon + ': %.2f\t%.3f' % (mean(freqList), std(freqList)))

            meanRbsMotifDic = {}
            for rbsMotifDic in data['rbs_motif']:
                for motif, freq in rbsMotifDic.iteritems():
                    meanRbsMotifDic[motif] = meanRbsMotifDic.get(motif, []) + [freq]

            for motif, freqList in meanRbsMotifDic.iteritems():
                fout.write('\t' + motif + ': %.2f\t%.3f' % (mean(freqList), std(freqList)))
            fout.write('\n')

        fout.close()

    def run(self):
        img = IMG()

        # get genome ids of all prokaryotes and euks
        print 'Reading Proks and Euks metadata.'
        bHeader = True
        genomeIdToGroup = {}
        missingGenomeData = {}
        for line in open(img.metadataFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            genomeId = lineSplit[0].strip()
            domain = lineSplit[1].strip()

            if os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + '.fna'):
                if domain == 'Bacteria' or domain == 'Archaea':
                    genomeIdToGroup[genomeId] = 'd__' + domain
                elif domain == 'Eukaryota':
                    phylum = lineSplit[6].strip()

                    if phylum == 'Apicomplexa' or phylum == 'Arthropoda' or phylum == 'Ascomycota' or phylum == 'Chlorophyta' or phylum == 'Chordata' or phylum == 'Streptophyta':
                        genomeIdToGroup[genomeId] = 'p__' + phylum
            else:
                missingGenomeData[domain] = missingGenomeData.get(domain, 0) + 1

        # get genome ids of all viruses
        print 'Reading Virus metadata.'
        bHeader = True
        for line in open(img.virusMetadataFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            genomeId = lineSplit[0].strip()
            domain = lineSplit[1].strip()

            if os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + '.fna'):
                genomeIdToGroup[genomeId] = 'd__' + domain
            else:
                missingGenomeData[domain] = missingGenomeData.get(domain, 0) + 1

        # report results
        print 'Number of valid genomes: ' + str(len(genomeIdToGroup))
        print 'Number of genomes missing genomic data: '
        for domain, count in missingGenomeData.iteritems():
            print '  ' + domain + ': ' + str(count)

        # process all genomes
        charByGroup = {}
        for genomeId, groupName in genomeIdToGroup.iteritems():
            print genomeId
            self.runProdigal(genomeId)
            self.getCharacteristics(genomeId, groupName, charByGroup)

        # write out results
        self.writeResults(charByGroup)

if __name__ == '__main__':
    genomeCharacteristics = GenomeCharacteristics()
    genomeCharacteristics.run()
