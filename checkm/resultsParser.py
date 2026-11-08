###############################################################################
#
# resultsParser.py - Parse and output results.
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

import sys
import os
import ast
import copy
from collections import defaultdict
import logging

import checkm.prettytable as prettytable

from checkm.defaultValues import DefaultValues
from checkm.common import reassignStdOut, restoreStdOut, checkFileExists
from checkm.coverage import Coverage
from checkm.hmmer import HMMERParser, FormatError, HMMERError

from checkm.util.pfam import PFAM
from checkm.util.seqUtils import readFasta


class ResultsParser():
    """Parse output of Prodigal+HMMER run and derived statistics."""

    def __init__(self, binIdToModels):
        self.logger = logging.getLogger('timestamp')

        self.results = {}
        self.models = binIdToModels

    def analyseResults(self,
                       outDir,
                       binStatsFile,
                       hmmTableFile,
                       bIgnoreThresholds=False,
                       evalueThreshold=DefaultValues.E_VAL,
                       lengthThreshold=DefaultValues.LENGTH,
                       bSkipPseudoGeneCorrection=False,
                       bSkipAdjCorrection=False,
                       ):
        """Parse results in the output directory."""

        # read bin and sequence stats into dictionaries
        binStats = self.parseBinStats(outDir, binStatsFile)

        # get hits to each bin
        self.parseBinHits(outDir, hmmTableFile, bSkipAdjCorrection, bIgnoreThresholds,
                          evalueThreshold, lengthThreshold, bSkipPseudoGeneCorrection, binStats)

        return binStats

    def cacheResults(self, outDir, binIdToBinMarkerSets, bIndividualMarkers):
        # cache critical results to file
        self._writeBinStatsExt(outDir, binIdToBinMarkerSets, bIndividualMarkers)
        self._writeMarkerGeneStats(outDir, binIdToBinMarkerSets, bIndividualMarkers)

    def parseBinHits(self, outDir,
                     hmmTableFile,
                     bSkipAdjCorrection=False,
                     bIgnoreThresholds=False,
                     evalueThreshold=DefaultValues.E_VAL,
                     lengthThreshold=DefaultValues.LENGTH,
                     bSkipPseudoGeneCorrection=False,
                     binStats=None):
        """ Parse HMM hits for each bin."""

        if not self.models:
            self.logger.error(
                'Models must be parsed before identifying HMM hits.')
            sys.exit(1)

        self.logger.info('Parsing HMM hits to marker genes:')

        numBinsProcessed = 0
        for binId in self.models:
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numBinsProcessed += 1
                statusStr = '    Finished parsing hits for %d of %d (%.2f%%) bins.' % (
                    numBinsProcessed, len(self.models), float(numBinsProcessed) * 100 / len(self.models))
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

            if binStats != None:
                resultsManager = ResultsManager(
                    binId, self.models[binId], bIgnoreThresholds, evalueThreshold, lengthThreshold, bSkipPseudoGeneCorrection, binStats[binId])
            else:
                resultsManager = ResultsManager(
                    binId, self.models[binId], bIgnoreThresholds, evalueThreshold, lengthThreshold, bSkipPseudoGeneCorrection)

            hmmerTableFile = os.path.join(outDir, 'bins', binId, hmmTableFile)
            self.parseHmmerResults(
                hmmerTableFile, resultsManager, bSkipAdjCorrection)
            self.results[binId] = resultsManager

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

    def _writeBinStatsExt(self, directory, binIdToBinMarkerSets, bIndividualMarkers):
        binStatsExtFile = os.path.join(directory, 'storage', DefaultValues.BIN_STATS_EXT_OUT)
        
        fout = open(binStatsExtFile, 'w')

        for binId in sorted(self.results):
            binStatsExt = self.results[binId].getSummary(
                binIdToBinMarkerSets[binId], bIndividualMarkers, outputFormat=2)
            binStatsExt.update(self.results[binId].geneCopyNumber(
                binIdToBinMarkerSets[binId]))
            fout.write(binId + '\t' + str(binStatsExt) + '\n')
        fout.close()

    def _writeMarkerGeneStats(self, directory, binIdToBinMarkerSets, bIndividualMarkers):
        markerGenesFile = os.path.join(directory, 'storage', DefaultValues.MARKER_GENE_STATS)

        fout = open(markerGenesFile, 'w')

        for binId in sorted(self.results):
            markerGenesHits = self.results[binId].getSummary(
                binIdToBinMarkerSets[binId], bIndividualMarkers, outputFormat=8)
            fout.write(binId + '\t' + str(markerGenesHits) + '\n')
        fout.close()

    def parseBinStats(self, resultsFolder, binStatsFile):
        """Read bin statistics from file."""
        binStatsFile = os.path.join(resultsFolder, 'storage', binStatsFile)

        checkFileExists(binStatsFile)

        binStats = {}
        with open(binStatsFile, 'r') as f:
            for line in f:
                lineSplit = line.split('\t')
                gid = lineSplit[0]
                data = lineSplit[1]
                binStats[gid] = ast.literal_eval(data)

        return binStats

    def parseBinStatsExt(self, resultsFolder):
        """Read bin statistics from file."""
        binStatsExtFile = os.path.join(
            resultsFolder, 'storage', DefaultValues.BIN_STATS_EXT_OUT)

        checkFileExists(binStatsExtFile)

        binStatsExt = {}
        with open(binStatsExtFile, 'r') as f:
            for line in f:
                lineSplit = line.split('\t')
                binStatsExt[lineSplit[0]] = ast.literal_eval(lineSplit[1])

        return binStatsExt

    def parseMarkerGeneStats(self, resultsFolder):
        """Read bin statistics from file."""
        markerGeneStatsFile = os.path.join(
            resultsFolder, 'storage', DefaultValues.MARKER_GENE_STATS)

        checkFileExists(markerGeneStatsFile)

        markerGeneStats = {}
        with open(markerGeneStatsFile, 'r') as f:
            for line in f:
                lineSplit = line.split('\t')
                markerGeneStats[lineSplit[0]] = ast.literal_eval(lineSplit[1])

        return markerGeneStats

    def parseHmmerResults(self, fileName, resultsManager, bSkipAdjCorrection):
        """Parse HMMER results."""
        try:
            with open(fileName, 'r') as hmmerHandle:
                HP = HMMERParser(hmmerHandle)

                while True:
                    hit = HP.next()
                    if hit is None:
                        break
                    resultsManager.addHit(hit)
        except IOError as detail:
            # a missing table must not be reported as a bin without marker genes
            self.logger.error('Failed to read HMMER results for bin %s: %s' % (resultsManager.binId, detail))
            sys.exit(1)
        except (Exception, FormatError, HMMERError) as detail:
            # (FormatError and HMMERError derive from BaseException)
            self.logger.error('Failed to parse HMMER results in %s: %s: %s' % (fileName, type(detail).__name__, detail))
            sys.exit(1)

        # retain only best hit to PFAM clan
        pfam = PFAM(DefaultValues.PFAM_CLAN_FILE)
        resultsManager.markerHits = pfam.filterHitsFromSameClan(
            resultsManager.markerHits)

        # correct for errors in ORF calling
        if not bSkipAdjCorrection:
            resultsManager.identifyAdjacentMarkerGenes()

    def _getHeader(self, outputFormat, binMarkerSets, coverageBinProfiles=None, table=None):
        """Get header for requested output table."""

        # keep count of single, double, triple genes etc...
        if outputFormat == 1:
            header = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', '0',
                      '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        elif outputFormat == 2:
            header = ['Bin Id']
            header += ['Marker lineage', '# genomes',
                       '# markers', '# marker sets']
            header += ['Completeness', 'Contamination', 'Strain heterogeneity']
            header += ['Genome size (bp)', '# ambiguous bases',
                       '# scaffolds', '# contigs']
            header += ['N50 (scaffolds)', 'N50 (contigs)',
                       'Mean scaffold length (bp)', 'Mean contig length (bp)']
            header += ['Longest scaffold (bp)', 'Longest contig (bp)']
            header += ['GC', 'GC std (scaffolds > 1kbp)']
            header += ['Coding density',
                       'Translation table', '# predicted genes']
            header += ['0', '1', '2', '3', '4', '5+']

            if coverageBinProfiles != None:
                for bamId in coverageBinProfiles[list(coverageBinProfiles.keys())[0]]:
                    header += ['Coverage (' + bamId + ')',
                               'Coverage std (' + bamId + ')']

            if DefaultValues.MIN_SEQ_LEN_GC_STD != 1000:
                self.logger.error('Labeling error: GC std (scaffolds > 1kbp)')
                sys.exit(1)
        elif outputFormat == 3:
            header = ['Bin Id', 'Node Id', 'Marker lineage', '# genomes', '# markers', '# marker sets',
                      '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        elif outputFormat == 4:
            header = None
        elif outputFormat == 5:
            header = ['Bin Id', 'Marker Id', 'Gene Id']
        elif outputFormat == 6:
            header = ['Bin Id', 'Marker Id', 'Gene Ids']
        elif outputFormat == 7:
            header = ['Bin Id', 'Marker Id', 'Gene Ids']
        elif outputFormat == 8:
            header = ['Bin Id', 'Gene Id',
                      '{Marker Id, Start position, End position}']
        elif outputFormat == 9:
            header = ['Bin Id', 'Contig', 'Gene Number', 'Gene Start', 'Gene End',
                      'Gene Strand', 'Prot Length', 'Marker Id', 'Align Start', 'Align End', 'Sequence']

        return header

    def printSummary(self, outputFormat, aai, binIdToBinMarkerSets, bIndividualMarkers, coverageFile, bTabTable, outFile, anaFolder):
        # redirect output
        oldStdOut = reassignStdOut(outFile)

        coverageBinProfiles = None
        if coverageFile:
            coverage = Coverage(1)
            coverageBinProfiles = coverage.binProfiles(coverageFile)

        prettyTableFormats = [1, 2, 3, 9]

        header = self._getHeader(outputFormat, binIdToBinMarkerSets[list(binIdToBinMarkerSets.keys())[0]], coverageBinProfiles, bTabTable)

        if bTabTable or outputFormat not in prettyTableFormats:
            bTabTable = True
            pTable = None

            if header != None:
                print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.2'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        seqsReported = 0
        for binId in sorted(self.results.keys()):
            seqsReported += self.results[binId].printSummary(
                outputFormat, aai, binIdToBinMarkerSets[binId], bIndividualMarkers, coverageBinProfiles, pTable, anaFolder)

        if outputFormat in [6, 7] and seqsReported == 0:
            print('[No marker genes satisfied the reporting criteria.]')

        if not bTabTable:
            if outputFormat in [1, 2]:
                print(pTable.get_string(sortby='Completeness', reversesort=True))
            else:
                # only print if there are rows
                if pTable.get_string(print_empty=False):
                    print(pTable.get_string(print_empty=False))

        # restore stdout
        restoreStdOut(outFile, oldStdOut)


class ResultsManager():
    """Store all the results for a single bin"""

    def __init__(self, binId, models,
                 bIgnoreThresholds=False,
                 evalueThreshold=DefaultValues.E_VAL,
                 lengthThreshold=DefaultValues.LENGTH,
                 bSkipPseudoGeneCorrection=False,
                 binStats=None):
        self.binId = binId
        self.markerHits = {}
        self.bIgnoreThresholds = bIgnoreThresholds
        self.evalueThreshold = evalueThreshold
        self.lengthThreshold = lengthThreshold
        self.bSkipPseudoGeneCorrection = bSkipPseudoGeneCorrection
        self.models = models
        self.binStats = binStats

    def vetHit(self, hit):
        """Check if hit meets required thresholds."""
        model = self.models[hit.query_accession]

        # preferentially use model specific bit score thresholds, before
        # using the user specified e-value and length criteria

        # Give preference to the gathering threshold unless the model
        # is marked as TIGR (i.e., TIGRFAM model)

        if not self.bSkipPseudoGeneCorrection:
            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length / float(hit.query_length)
            if length_perc < DefaultValues.PSEUDOGENE_LENGTH:
                return False

        if model.nc != None and not self.bIgnoreThresholds and 'TIGR' in model.acc:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        elif model.ga != None and not self.bIgnoreThresholds:
            if model.ga[0] <= hit.full_score and model.ga[1] <= hit.dom_score:
                return True
        elif model.tc != None and not self.bIgnoreThresholds:
            if model.tc[0] <= hit.full_score and model.tc[1] <= hit.dom_score:
                return True
        elif model.nc != None and not self.bIgnoreThresholds:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        else:
            if hit.full_e_value > self.evalueThreshold:
                return False

            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length / float(hit.query_length)
            if length_perc >= self.lengthThreshold:
                return True

        return False

    def addHit(self, hit):
        """Process hit and add it to the set of markers if it passes filtering criteria."""
        if self.vetHit(hit):
            if hit.query_accession in self.markerHits:
                # retain only the best domain hit for a given marker to a specific ORF
                previousHitToORF = None
                for h in self.markerHits[hit.query_accession]:
                    if h.target_name == hit.target_name:
                        previousHitToORF = h
                        break

                if not previousHitToORF:
                    self.markerHits[hit.query_accession].append(hit)
                else:
                    if previousHitToORF.dom_score < hit.dom_score:
                        self.markerHits[hit.query_accession].append(hit)
                        self.markerHits[hit.query_accession].remove(
                            previousHitToORF)

            else:
                self.markerHits[hit.query_accession] = [hit]

    def identifyAdjacentMarkerGenes(self):
        """Identify adjacent marker genes and exclude these from the contamination estimate."""

        # check for adjacent ORFs with hits to the same marker gene
        for markerId, hits in self.markerHits.items():

            bCombined = True
            while bCombined:
                bCombined = False
                for i in range(0, len(hits)):
                    orfI = hits[i].target_name
                    scaffoldIdI = orfI[0:orfI.rfind('_')]

                    for j in range(i + 1, len(hits)):
                        orfJ = hits[j].target_name
                        scaffoldIdJ = orfJ[0:orfJ.rfind('_')]

                        # check if hits are on adjacent ORFs
                        if scaffoldIdI == scaffoldIdJ:
                            try:
                                orfNumI = int(orfI[orfI.rfind('_') + 1:])
                                orfNumJ = int(orfJ[orfJ.rfind('_') + 1:])
                            except:
                                # it appears called genes are not labeled
                                # according to the prodigal format, so
                                # it is not possible to perform this correction
                                break

                            if abs(orfNumI - orfNumJ) == 1:
                                # Hits to adjacent ORFs are combined. If the hits are to
                                # different parts of the HMM, the marker gene was most likely
                                # split by a gene calling or assembly error. If the hits
                                # overlap on the HMM this could represent contamination,
                                # but it seems more likely that adjacent genes hitting
                                # the same marker represent legitimate gene duplication,
                                # a gene calling error, or an assembly error and thus
                                # should not be treated as contamination.
                                bCombined = True
                                break

                    if bCombined:
                        # copy the hit so the original hit to each ORF is retained
                        # (used to report the position of each gene, e.g. qa -o 9)
                        newHit = copy.copy(hits[i])
                        newHit.component_hits = (getattr(hits[i], 'component_hits', [hits[i]])
                                                 + getattr(hits[j], 'component_hits', [hits[j]]))

                        # produce concatenated label indicating the two genes being combined
                        orfA, orfB = sorted([orfI, orfJ])
                        newHit.target_name = DefaultValues.SEQ_CONCAT_CHAR.join([
                                                                                orfA, orfB])

                        newHit.target_length = hits[i].target_length + \
                            hits[j].target_length
                        newHit.hmm_from = min(
                            hits[i].hmm_from, hits[j].hmm_from)
                        newHit.hmm_to = max(hits[i].hmm_to, hits[j].hmm_to)

                        newHit.ali_from = min(
                            hits[i].ali_from, hits[j].ali_from)
                        newHit.ali_to = max(hits[i].ali_to, hits[j].ali_to)

                        newHit.env_from = min(
                            hits[i].env_from, hits[j].env_from)
                        newHit.env_to = max(hits[i].env_to, hits[j].env_to)

                        hits.remove(hits[j])
                        hits.remove(hits[i])

                        hits.append(newHit)

                        break

            self.markerHits[markerId] = hits

    def countUniqueHits(self):
        """Determine number of unique and multiple hits."""
        uniqueHits = 0
        multiCopyHits = 0
        for hits in self.markerHits.values():
            if len(hits) == 1:
                uniqueHits += 1
            elif len(hits) > 1:
                multiCopyHits += 1

        return uniqueHits, multiCopyHits

    def hitsToMarkerGene(self, markerSet):
        """Determine number of times each marker gene is hit."""
        ret = dict()
        for marker in sorted(markerSet.getMarkerGenes()):
            # use get() so absent markers are not inserted into the hit dictionary
            ret[marker] = len(self.markerHits.get(marker, []))

        return ret

    def geneCountsForSelectedMarkerSet(self, binMarkerSets, bIndividualMarkers):
        """Report number of times each marker gene was found along with
           a completeness and contamination estimation for the selected marker
           set associated with a bin."""
        geneCounts = self.geneCounts(
            binMarkerSets.selectedMarkerSet(), self.markerHits, bIndividualMarkers)

        return geneCounts

    def geneCounts(self, markerSet, markerHits, bIndividualMarkers):
        """ Determine number of marker genes with 0-5 hits
            as well as the total completeness and contamination."""
        geneCounts = [0] * 6
        for marker in markerSet.getMarkerGenes():
            if marker in markerHits:
                if len(markerHits[marker]) > 5:
                    markerCount = 5
                else:
                    markerCount = len(markerHits[marker])
            else:
                markerCount = 0

            geneCounts[markerCount] += 1

        percComp, percCont = markerSet.genomeCheck(
            markerHits, bIndividualMarkers)

        geneCounts.append(percComp)
        geneCounts.append(percCont)

        return geneCounts

    def geneCopyNumber(self, binMarkerSets):
        """ Determine number of times each marker gene is present."""
        geneCopyNumber = {}
        geneCopyNumber['GCN0'] = []
        geneCopyNumber['GCN1'] = []
        geneCopyNumber['GCN2'] = []
        geneCopyNumber['GCN3'] = []
        geneCopyNumber['GCN4'] = []
        geneCopyNumber['GCN5+'] = []

        markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

        for marker in sorted(self.models):
            if marker not in markerGenes:
                continue

            markerId = os.path.splitext(marker)[0]
            if marker in self.markerHits:
                if len(self.markerHits[marker]) >= 5:
                    geneCopyNumber['GCN5+'].append(markerId)
                else:
                    geneCopyNumber['GCN' +
                                   str(len(self.markerHits[marker]))].append(markerId)
            else:
                geneCopyNumber['GCN0'].append(markerId)

        return geneCopyNumber

    def getSummary(self, binMarkerSets, bIndividualMarkers, outputFormat):
        """Get dictionary containing information about bin.

        Output format 2 gives the extended bin statistics (bin_stats_ext.tsv) and
        output format 8 the position of marker genes (marker_gene_stats.tsv).
        """
        summary = {}

        if outputFormat == 2:
            selectedMarkerSet = binMarkerSets.selectedMarkerSet()
            data = self.geneCountsForSelectedMarkerSet(
                binMarkerSets, bIndividualMarkers)

            summary['marker lineage'] = selectedMarkerSet.lineageStr
            summary['# genomes'] = selectedMarkerSet.numGenomes
            summary['# markers'] = selectedMarkerSet.numMarkers()
            summary['# marker sets'] = selectedMarkerSet.numSets()
            summary['0'] = data[0]
            summary['1'] = data[1]
            summary['2'] = data[2]
            summary['3'] = data[3]
            summary['4'] = data[4]
            summary['5+'] = data[5]
            summary['Completeness'] = data[6]
            summary['Contamination'] = data[7]
            summary.update(self.binStats)

        elif outputFormat == 8:
            # tabular - print only position of marker genes
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

            genesWithMarkers = {}
            for marker, hit_list in sorted(self.markerHits.items()):
                if marker not in markerGenes:
                    continue

                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(
                        hit.target_name, []) + [hit]

            for geneId, hits in sorted(genesWithMarkers.items()):
                summary[geneId] = {}
                for hit in hits:
                    summary[geneId][hit.query_accession] = summary[geneId].get(
                        hit.query_accession, []) + [[hit.ali_from, hit.ali_to]]

        return summary

    def printSummary(self, outputFormat, aai, binMarkerSets, bIndividualMarkers, coverageBinProfiles=None, table=None, anaFolder=None):
        """Print out information about bin."""
        if outputFormat == 1:
            selectedMarkerSet = binMarkerSets.selectedMarkerSet()

            lineageStr = selectedMarkerSet.lineageStr
            if str(selectedMarkerSet.UID) != '0':
                # marker sets from a file of HMMs have the UID 0
                lineageStr += ' (' + str(selectedMarkerSet.UID) + ')'

            data = self.geneCountsForSelectedMarkerSet(
                binMarkerSets, bIndividualMarkers)
            row = "%s\t%s\t%d\t%d\t%d\t%s\t%0.2f\t%0.2f\t%0.2f" % (self.binId, lineageStr,
                                                                   selectedMarkerSet.numGenomes, selectedMarkerSet.numMarkers(), selectedMarkerSet.numSets(),
                                                                   "\t".join(
                                                                       [str(data[i]) for i in range(6)]),
                                                                   data[6],
                                                                   data[7],
                                                                   aai.aaiMeanBinHetero.get(
                                                                       self.binId, 0.0)
                                                                   )
            if table == None:
                print(row)
            else:
                table.add_row([self.binId, lineageStr, selectedMarkerSet.numGenomes, selectedMarkerSet.numMarkers(
                ), selectedMarkerSet.numSets()] + data + [aai.aaiMeanBinHetero.get(self.binId, 0.0)])
        elif outputFormat == 2:
            selectedMarkerSet = binMarkerSets.selectedMarkerSet()

            lineageStr = selectedMarkerSet.lineageStr
            if str(selectedMarkerSet.UID) != '0':
                # marker sets from a file of HMMs have the UID 0
                lineageStr += ' (' + str(selectedMarkerSet.UID) + ')'

            data = self.geneCountsForSelectedMarkerSet(
                binMarkerSets, bIndividualMarkers)

            if table == None:
                row = self.binId
                row += '\t%s\t%d\t%d\t%d' % (lineageStr, selectedMarkerSet.numGenomes,
                                             selectedMarkerSet.numMarkers(), selectedMarkerSet.numSets())
                row += '\t%0.2f\t%0.2f\t%0.2f' % (
                    data[6], data[7], aai.aaiMeanBinHetero.get(self.binId, 0.0))
                row += '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (self.binStats['Genome size'], self.binStats['# ambiguous bases'],
                                                                     self.binStats['# scaffolds'], self.binStats['# contigs'],
                                                                     self.binStats['N50 (scaffolds)'], self.binStats['N50 (contigs)'],
                                                                     self.binStats['Mean scaffold length'], self.binStats['Mean contig length'],
                                                                     self.binStats['Longest scaffold'], self.binStats['Longest contig'])
                row += '\t%.1f\t%.2f' % (self.binStats['GC']
                                         * 100, self.binStats['GC std'] * 100)
                row += '\t%.2f\t%d\t%d' % (self.binStats['Coding density'] * 100,
                                           self.binStats['Translation table'], self.binStats['# predicted genes'])
                row += '\t' + '\t'.join([str(data[i]) for i in range(6)])

                if coverageBinProfiles:
                    if self.binId in coverageBinProfiles:
                        for _, coverageStats in coverageBinProfiles[self.binId].items():
                            row += '\t%.2f\t%.2f' % (
                                coverageStats[0], coverageStats[1])
                    else:
                        for bamId in coverageBinProfiles[list(coverageBinProfiles.keys())[0]]:
                            row += '\t%.2f\t%.2f' % (0, 0)
                print(row)
            else:
                row = [self.binId, lineageStr, selectedMarkerSet.numGenomes,
                       selectedMarkerSet.numMarkers(), selectedMarkerSet.numSets()]
                row.extend(
                    [data[6], data[7], aai.aaiMeanBinHetero.get(self.binId, 0.0)])
                row.extend([self.binStats['Genome size'], self.binStats['# ambiguous bases'], self.binStats['# scaffolds'],
                            self.binStats['# contigs'], self.binStats['N50 (scaffolds)'], self.binStats[
                    'N50 (contigs)'],
                    int(self.binStats['Mean scaffold length']), int(
                                self.binStats['Mean contig length']),
                    self.binStats['Longest scaffold'], self.binStats['Longest contig']])
                row.extend([self.binStats['GC'] * 100,
                           self.binStats['GC std'] * 100])
                row.extend([self.binStats['Coding density'] * 100,
                           self.binStats['Translation table'], self.binStats['# predicted genes']])
                row.extend(data[0:6])

                if coverageBinProfiles:
                    if self.binId in coverageBinProfiles:
                        for _, coverageStats in coverageBinProfiles[self.binId].items():
                            row.extend(coverageStats)
                    else:
                        for bamId in coverageBinProfiles[list(coverageBinProfiles.keys())[0]]:
                            row.extend([0, 0])

                table.add_row(row)
        elif outputFormat == 3:
            for ms in binMarkerSets.markerSetIter():
                data = self.geneCounts(ms, self.markerHits, bIndividualMarkers)
                row = "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%0.2f\t%0.2f\t%0.2f" % (self.binId, ms.UID, ms.lineageStr, ms.numGenomes,
                                                                           ms.numMarkers(), ms.numSets(),
                                                                           "\t".join(
                                                                               [str(data[i]) for i in range(6)]),
                                                                           data[6],
                                                                           data[7],
                                                                           aai.aaiMeanBinHetero.get(
                                                                               self.binId, 0.0)
                                                                           )
                if table == None:
                    print(row)
                else:
                    table.add_row([self.binId, ms.UID, ms.lineageStr, ms.numGenomes, ms.numMarkers(
                    ), ms.numSets()] + data + [aai.aaiMeanBinHetero.get(self.binId, 0.0)])

        elif outputFormat == 4:
            selectedMarkerSet = binMarkerSets.selectedMarkerSet()
            data = self.hitsToMarkerGene(binMarkerSets.selectedMarkerSet())
            row = "Node Id: %s; Marker lineage: %s" % (
                selectedMarkerSet.UID, selectedMarkerSet.lineageStr)
            for marker in data:
                row += '\t' + marker
            print(row)

            row = self.binId
            for count in data.values():
                row += '\t' + str(count)
            print(row)

            print()
        elif outputFormat == 5:
            # tabular of bin_id, marker, contig_id
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

            for marker, hit_list in sorted(self.markerHits.items()):
                if marker not in markerGenes:
                    continue

                for geneId in sorted(hit.target_name for hit in hit_list):
                    print(self.binId, marker, geneId, sep='\t', end='\n')

        elif outputFormat == 6:
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

            seqsReported = 0
            for marker, hitList in sorted(self.markerHits.items()):
                if marker not in markerGenes:
                    continue

                if len(hitList) >= 2:
                    print(self.binId, marker, sep='\t', end='\t')

                    scaffoldIds = []
                    for hit in hitList:
                        scaffoldIds.append(hit.target_name)

                    print(','.join(sorted(scaffoldIds)), end='\n')

                    seqsReported += 1

            return seqsReported

        elif outputFormat == 7:
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

            seqsReported = 0
            for marker, hitList in sorted(self.markerHits.items()):
                if marker not in markerGenes:
                    continue

                if len(hitList) >= 2:
                    scaffoldsWithMultipleHits = set()
                    for i in range(0, len(hitList)):
                        scaffoldId = hitList[i].target_name[0:hitList[i].target_name.rfind(
                            '_')]
                        for j in range(i + 1, len(hitList)):
                            if scaffoldId == hitList[j].target_name[0:hitList[j].target_name.rfind('_')]:
                                scaffoldsWithMultipleHits.add(
                                    hitList[i].target_name)
                                scaffoldsWithMultipleHits.add(
                                    hitList[j].target_name)

                    if len(scaffoldsWithMultipleHits) >= 2:
                        print(self.binId, marker, sep='\t', end='\t')
                        print(
                            ','.join(sorted(list(scaffoldsWithMultipleHits))), end='\n')
                        seqsReported += 1

            return seqsReported

        elif outputFormat == 8:
            # tabular - print only position of marker genes
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

            genesWithMarkers = {}
            for marker, hit_list in sorted(self.markerHits.items()):
                if marker not in markerGenes:
                    continue

                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(
                        hit.target_name, []) + [hit]

            for geneId, hits in sorted(genesWithMarkers.items()):
                rowStr = self.binId + '\t' + geneId
                for hit in hits:
                    rowStr += '\t' + hit.query_accession + ',' + \
                        str(hit.ali_from) + ',' + str(hit.ali_to)
                print(rowStr)

        # Hunter Cameron, May 29, 2015 - print a fasta of marker genes
        elif outputFormat == 9:
            # tabular of bin_id, marker, contig_id

            # check for the analyze folder for later use
            if anaFolder is None:
                raise ValueError(
                    "AnaFolder must not be None for outputFormat 9")

            # ## build a dict to link genes with marker gene alignment information
            # (a gene may be hit by several markers, and hits to adjacent genes
            # may have been combined into a single hit, e.g. 'geneA&&geneB')
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()
            hitInfo = defaultdict(list)
            for marker in sorted(self.markerHits):
                if marker not in markerGenes:
                    continue

                for hit in self.markerHits[marker]:
                    for geneHit in getattr(hit, 'component_hits', [hit]):
                        for geneId in geneHit.target_name.split(DefaultValues.SEQ_CONCAT_CHAR):
                            hitInfo[geneId].append({
                                "marker": marker,
                                "ali_from": str(geneHit.ali_from),
                                "ali_to": str(geneHit.ali_to)
                            })

            # ## Open genes.faa and print the ones that were found with some descriptive info in the header
            path_to_genes = os.path.join(anaFolder, "bins", self.binId, DefaultValues.PRODIGAL_AA)

            # get only the seqs we need and their information as a dict
            seqs = readFasta(path_to_genes, trimHeader=False)

            def parse_header(header):
                """Parse gene information from a Prodigal header.

                Headers of user-supplied protein files (--genes) may not follow the
                Prodigal format, in which case 'NA' is reported for unknown fields.
                """
                elems = header.split(" # ")
                fields = elems[0].split(None, 1)
                gene_name = fields[0] if fields else ''

                # remove the gene number from Prodigal to get the original contig name
                contig_name = gene_name
                gene_num = 'NA'
                if '_' in gene_name:
                    ctg_name, num = gene_name.rsplit("_", 1)
                    if num.isdigit():
                        contig_name, gene_num = ctg_name, num

                if len(elems) >= 4:
                    gene_start, gene_end, gene_strand = [e.strip() for e in elems[1:4]]
                else:
                    gene_start, gene_end, gene_strand = 'NA', 'NA', 'NA'

                return gene_name, contig_name, gene_num, gene_start, gene_end, gene_strand

            genes = []
            for header in seqs.keys():
                gene_info = parse_header(header)
                if gene_info[0] in hitInfo:
                    genes.append((header, gene_info))

            def sort_genes(gene):
                """ sorts genes by contig and gene number """
                gene_name, contig_name, gene_num = gene[1][0:3]
                if gene_num.isdigit():
                    return contig_name, int(gene_num), gene_name
                return contig_name, -1, gene_name

            for header, (gene_name, contig_name, gene_num, gene_start, gene_end, gene_strand) in sorted(genes, key=sort_genes):
                seq = seqs[header]
                for marker_hit in hitInfo[gene_name]:
                    # if table output not specified, print FASTA
                    if table != None:
                        gene_info = "geneId={};start={};end={};strand={};protlen={}".format(
                            gene_num, gene_start, gene_end, gene_strand, str(len(seq)))

                        marker_info = "marker={};mstart={};mend={}".format(
                            marker_hit["marker"],
                            marker_hit["ali_from"],
                            marker_hit["ali_to"])

                        # new header will be the bin name, contig name, gene info, and marker info separated by spaces
                        new_header = ">" + \
                            " ".join([self.binId, contig_name,
                                     gene_info, marker_info])

                        print(new_header, seq, sep="\n")
                    # otherwise, print a table
                    else:
                        print("\t".join([
                            self.binId,
                            contig_name,
                            gene_num,
                            gene_start,
                            gene_end,
                            gene_strand,
                            str(len(seq)),
                            marker_hit["marker"],
                            marker_hit["ali_from"],
                            marker_hit["ali_to"],
                            seq
                        ]))

        return 0
