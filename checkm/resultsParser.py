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

from __future__ import print_function

import sys
import os 
import ast
from collections import defaultdict
import logging

import prettytable

import defaultValues 
from common import reassignStdOut, restoreStdOut, checkFileExists
from coverage import Coverage

from hmmer import HMMERParser
from hmmerModelParser import HmmModelParser

from lib.pfam import PFAM

class ResultsParser():
    """Parse output of Prodigal+HMMER run and derived statistics."""
    def __init__(self):
        self.logger = logging.getLogger()

        self.results = {}
        self.models = defaultdict(dict)
    
    def analyseResults(self,
                       outDir,
                       binStatsFile,
                       seqStatsFile,
                       hmmTableFile,
                       binIdToHmmModelFile,
                       bIgnoreThresholds = False,
                       evalueThreshold = defaultValues.E_VAL,
                       lengthThreshold = defaultValues.LENGTH,
                       bSkipOrfCorrection = False,
                       ):
        """Parse results in the output directory."""

        # parse information from HMMs
        self.parseHmmerModels(binIdToHmmModelFile)
        
        # read bin and sequence stats into dictionaries
        binStats = self.parseBinStats(outDir, binStatsFile)
        seqStats = self.parseSeqStats(outDir, seqStatsFile)

        # get hits to each bin
        self.parseBinHits(outDir, hmmTableFile, bSkipOrfCorrection, bIgnoreThresholds, evalueThreshold, lengthThreshold, binStats, seqStats)
            
    def cacheResults(self, outDir, binIdToBinMarkerSets, bIndividualMarkers):
        # cache critical results to file
        self.__writeBinStatsExt(outDir, binIdToBinMarkerSets, bIndividualMarkers)
        self.__writeMarkerGeneStats(outDir, binIdToBinMarkerSets, bIndividualMarkers)
        
    def parseHmmerModels(self, binIdToHmmModelFile):
        """Parse each bin's HMM file to collect model information."""
        for binId, hmmModelFile in binIdToHmmModelFile.iteritems():
            modelParser = HmmModelParser(hmmModelFile)
            for model in modelParser.parse():
                self.models[binId][model.acc] = model

    def parseBinHits(self, outDir, 
                     hmmTableFile, 
                     bSkipOrfCorrection = False, 
                     bIgnoreThresholds = False, 
                     evalueThreshold = defaultValues.E_VAL, 
                     lengthThreshold = defaultValues.LENGTH, 
                     binStats = None, 
                     seqStats = None):
        """ Parse HMM hits for each bin."""
        if not self.models:
            self.logger.error('  [Error] Models must be parsed before identifying HMM hits.')
            raise
        
        for binId in self.models:    
            if binStats != None and seqStats != None:
                resultsManager = ResultsManager(binId, self.models[binId], bIgnoreThresholds, evalueThreshold, lengthThreshold, binStats[binId], seqStats[binId])
            elif binStats == None and seqStats == None:
                resultsManager = ResultsManager(binId, self.models[binId], bIgnoreThresholds, evalueThreshold, lengthThreshold)
            else:
                self.logger.error('  [Error] Invalid parameter settings for binStats and seqStats.')
                raise
                   
            hmmerTableFile = os.path.join(outDir, 'bins', binId, hmmTableFile)
            self.parseHmmerResults(hmmerTableFile, resultsManager, bSkipOrfCorrection)
            self.results[binId] = resultsManager

    def __writeBinStatsExt(self, directory, binIdToBinMarkerSets, bIndividualMarkers):
        binStatsExt = {}
        for binId in self.results:
            binStatsExt[binId] = self.results[binId].getSummary(binIdToBinMarkerSets[binId], bIndividualMarkers, outputFormat=2)
            binStatsExt[binId].update(self.results[binId].geneCopyNumber())
           
        binStatsExtFile = os.path.join(directory, 'storage', defaultValues.BIN_STATS_EXT_OUT) 
        fout = open(binStatsExtFile, 'w')
        fout.write(str(binStatsExt))
        fout.close
        
    def __writeMarkerGeneStats(self, directory, binIdToBinMarkerSets, bIndividualMarkers):
        markerGenes = {}
        for binId in self.results:
            markerGenes[binId] = self.results[binId].getSummary(binIdToBinMarkerSets[binId], bIndividualMarkers, outputFormat=8)
             
        markerGenesFile = os.path.join(directory, 'storage', defaultValues.MARKER_GENE_STATS)
        fout = open(markerGenesFile, 'w')
        fout.write(str(markerGenes))
        fout.close
                         
    def parseBinStats(self, resultsFolder, binStatsFile):
        """Read bin statistics from file."""
        binStatsFile = os.path.join(resultsFolder, 'storage', binStatsFile)
        
        checkFileExists(binStatsFile)
        
        with open(binStatsFile, 'r') as f:
            s = f.read()
            binStats = ast.literal_eval(s)
            
        return binStats
    
    def parseBinStatsExt(self, resultsFolder):
        """Read bin statistics from file."""
        binStatsExtFile = os.path.join(resultsFolder, 'storage', defaultValues.BIN_STATS_EXT_OUT)
        
        checkFileExists(binStatsExtFile)
        
        with open(binStatsExtFile, 'r') as f:
            s = f.read()
            binStatsExt = ast.literal_eval(s)
            
        return binStatsExt
    
    def parseMarkerGeneStats(self, resultsFolder):
        """Read bin statistics from file."""
        markerGeneStatsFile = os.path.join(resultsFolder, 'storage', defaultValues.MARKER_GENE_STATS)
        
        checkFileExists(markerGeneStatsFile)
        
        with open(markerGeneStatsFile, 'r') as f:
            s = f.read()
            markerGeneStats = ast.literal_eval(s)
            
        return markerGeneStats
            
    def parseSeqStats(self, resultsFolder, seqStatsFile):
        """Read sequence statistics from file."""
        seqStatsFile = os.path.join(resultsFolder, 'storage', seqStatsFile)
        
        checkFileExists(seqStatsFile)
        
        with open(seqStatsFile, 'r') as f:
            s = f.read()
            seqStats = ast.literal_eval(s)
            
        return seqStats
            
    def parseHmmerResults(self, fileName, resultsManager, bSkipOrfCorrection):
        """Parse HMMER results."""
        try:
            with open(fileName, 'r') as hmmerHandle:
                try:
                    HP = HMMERParser(hmmerHandle)
                except:
                    print("Error opening HMM file: ", fileName)
                    raise
                
                while True:
                    hit = HP.next()
                    if hit is None:
                        break
                    resultsManager.addHit(hit)

            # retain only best hit to PFAM clan
            pfam = PFAM()
            resultsManager.markerHits = pfam.filterHitsFromSameClan(resultsManager.markerHits)
            
            # correct for errors in ORF calling      
            if not bSkipOrfCorrection:
                resultsManager.identifyOrfErrors()
                      
        except IOError as detail:
            sys.stderr.write(str(detail)+"\n")

    def __getHeader(self, outputFormat, coverageBinProfiles = None):
        """Get header for requested output table."""
                
        # keep count of single, double, triple genes etc...
        if outputFormat == 1:
            header = ['Bin Id','Marker lineage', '# markers','# marker sets','0','1','2','3','4','5+','Completeness','Contamination','Strain heterogeneity']
        elif outputFormat == 2:
            header = ['Bin Id']
            header += ['Marker lineage', '# markers', '# marker sets']
            header += ['Completeness','Contamination', 'Strain heterogeneity']
            header += ['Genome size (bp)', '# scaffolds', '# contigs', 'N50 (scaffolds)', 'N50 (contigs)', 'Longest scaffold (bp)', 'Longest contig (bp)']
            header += ['GC', 'GC std (scaffolds > 1Kbps)']
            header += ['Coding density (translation table 11)', '# predicted ORFs']
            header += ['0','1','2','3','4','5+']
            
            if coverageBinProfiles != None:
                for bamId in coverageBinProfiles[coverageBinProfiles.keys()[0]]:
                    header += ['Coverage (' + bamId + ')', 'Coverage std']

            if defaultValues.MIN_SEQ_LEN_GC_STD != 1000:
                self.logger.error('  [Error] Labeling error: GC std (scaffolds > 1Kbps)')
                sys.exit() 
        elif outputFormat == 3:
            header = ['Bin Id', 'Marker Id', 'Marker lineage', '# genomes', '# markers','# marker sets','0','1','2','3','4','5+','Completeness','Contamination','Strain heterogeneity']
        elif outputFormat == 4:
            header = []
        elif outputFormat == 5:
            header = [''] + self.models.keys()
        elif outputFormat == 6:
            header = ['Bin Id','Marker Id','Scaffold Id']
        elif outputFormat == 7:
            header = ['Bin Id','Marker Id','Scaffold Ids']
        elif outputFormat == 8:
            header = ['Bin Id','Scaffold Id','{Marker Id, Count}']
        elif outputFormat == 9:
            header = ['Bin Id','Gene Id','{Marker Id, Start position, End position}']
        elif outputFormat == 10:
            header = ['Scaffold Id','Bin Id','Length','# contigs','GC','# ORFs','Coding density','Marker Ids']
            
        return header
        
    def printSummary(self, outputFormat, aai, binIdToBinMarkerSets, bIndividualMarkers, coverageFile, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)
        
        coverageBinProfiles = None
        if coverageFile:
            coverage = Coverage(1)
            coverageProfile = coverage.parseCoverage(coverageFile)
            coverageBinProfiles = coverage.binProfiles(coverageProfile)
                
        prettyTableFormats = [1, 2, 3]      
          
        header = self.__getHeader(outputFormat, coverageBinProfiles) 
        if bTabTable or outputFormat not in prettyTableFormats: 
            bTabTable = True
            pTable = None
            print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.2'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        for binId in sorted(self.results.keys()):
            self.results[binId].printSummary(outputFormat, aai, binIdToBinMarkerSets[binId], bIndividualMarkers, coverageBinProfiles, pTable)
            
        if not bTabTable :  
            if outputFormat in [1,2]:
                print(pTable.get_string(sortby='Completeness', reversesort=True))
            else:
                print(pTable.get_string())
            
        # restore stdout   
        restoreStdOut(outFile, oldStdOut)     

class ResultsManager():
    """Store all the results for a single bin"""
    def __init__(self, binId, models,
                 bIgnoreThresholds = False, 
                 evalueThreshold = defaultValues.E_VAL, 
                 lengthThreshold = defaultValues.LENGTH, 
                 binStats=None, 
                 scaffoldStats=None):
        self.binId = binId
        self.markerHits = {}
        self.bIgnoreThresholds = bIgnoreThresholds
        self.evalueThreshold = evalueThreshold
        self.lengthThreshold = lengthThreshold
        self.models = models
        self.binStats = binStats
        self.scaffoldStats = scaffoldStats
    
    def vetHit(self, hit):
        """Check if hit meets required thresholds."""
        model = self.models[hit.query_accession]
        
        # preferentially use model specific bit score thresholds, before
        # using the user specified e-value and length criteria
        if model.isGatheringThreshold and not self.bIgnoreThresholds:
            if model.ga[0] <= hit.full_score and model.ga[1] <= hit.dom_score:
                return True
        elif model.isTrustedCutoff and not self.bIgnoreThresholds:
            if model.tc[0] <= hit.full_score and model.tc[1] <= hit.dom_score:
                return True
        elif model.isNoiseCutoff and not self.bIgnoreThresholds:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        else:
            if hit.full_e_value > self.evalueThreshold:
                return False

            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length/float(hit.query_length)
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
                        self.markerHits[hit.query_accession].remove(previousHitToORF)
                    
            else:
                self.markerHits[hit.query_accession] = [hit]
                
    def identifyOrfErrors(self):
        """Identify ORF errors affecting marker genes."""
 
        # check for adjacent ORFs with hits to the same marker gene
        for markerId, hits in self.markerHits.iteritems():
            
            bCombined = True
            while bCombined:
                for i in xrange(0, len(hits)):
                    orfI = hits[i].target_name
                    scaffoldIdI = orfI[0:orfI.rfind('_')]
                    
                    bCombined = False
                    for j in xrange(i+1, len(hits)):
                        orfJ = hits[j].target_name
                        scaffoldIdJ = orfJ[0:orfJ.rfind('_')]
                        
                        # check if hits are on adjacent ORFs
                        if scaffoldIdI == scaffoldIdJ:
                            orfNumI = int(orfI[orfI.rfind('_')+1:])
                            orfNumJ = int(orfJ[orfJ.rfind('_')+1:])
                            if abs(orfNumI - orfNumJ) == 1:
                                # check if hits are to different parts of the HMM
                                sI = hits[i].hmm_from
                                eI = hits[i].hmm_to
                                
                                sJ = hits[j].hmm_from
                                eJ = hits[j].hmm_to
    
                                if (sI <= sJ and eI > sJ) or (sJ <= sI and eJ > sI):
                                    # models overlap so treat these as unique hits
                                    # which may represent an assembly error or a true
                                    # gene duplication
                                    pass
                                else:
                                    # combine the two hits
                                    bCombined = True
                                    break
                                    
                    if bCombined:
                        newHit = hits[i]
                        newHit.target_name = defaultValues.SEQ_CONCAT_CHAR.join(sorted([orfI, orfJ]))
                        
                        newHit.target_length = hits[i].target_length + hits[j].target_length
                        
                        newHit.hmm_from = min(hits[i].hmm_from, hits[j].hmm_from)
                        newHit.hmm_to = min(hits[i].hmm_to, hits[j].hmm_to)
                        
                        newHit.ali_from = min(hits[i].ali_from, hits[j].ali_from)
                        newHit.ali_to = min(hits[i].ali_to, hits[j].ali_to)
                        
                        newHit.env_from = min(hits[i].env_from, hits[j].env_from)
                        newHit.env_to = min(hits[i].env_to, hits[j].env_to)
                        
                        hits.remove(hits[i])
                        hits.remove(hits[i])
                        
                        hits.append(newHit)

                        break
                    
            self.markerHits[markerId] = hits
                
    def hitsToMarkerGene(self, markerSet):
        """Determine number of times each marker gene is hit."""
        ret = dict()
        for marker in markerSet.getMarkerGenes():
            try:
                ret[marker] = len(self.markerHits[marker])
            except KeyError:
                ret[marker] = 0
                
        return ret

    def geneCountsForBestMarkerSet(self, binMarkerSets, bIndividualMarkers):
        """Report number of times each marker gene was found along with 
           a completeness and contamination estimation for the best marker
           set associated with a bin."""   
        bestPercComp = -1
        bestMarkerSet = None
        for ms in binMarkerSets.markerSetIter():
            geneCounts = self.geneCounts(ms, self.markerHits, bIndividualMarkers)
            
            if geneCounts[6] > bestPercComp:
                bestPercComp = geneCounts[6]
                bestGeneCount = geneCounts
                bestMarkerSet = ms
        
        return bestGeneCount, bestMarkerSet
        
    def geneCounts(self, markerSet, markerHits, bIndividualMarkers):
        """ Determine number of marker genes with 0-5 hits 
            as well as the total completeness and contamination."""   
        geneCounts = [0]*6
        multiCopyCount = 0
        for marker in markerSet.getMarkerGenes():
            if marker in markerHits:
                if len(markerHits[marker]) > 5:
                    markerCount = 5
                else:
                    markerCount = len(markerHits[marker])
                    
                multiCopyCount += (len(markerHits[marker]) - 1)
            else:
                markerCount = 0
            
            geneCounts[markerCount] += 1
            
        percComp, percCont = markerSet.genomeCheck(markerHits, bIndividualMarkers)
        
        geneCounts.append(percComp)
        geneCounts.append(percCont)
        
        return geneCounts
    
    def geneCopyNumber(self):
        """ Determine number of times each marker gene is present."""
        geneCopyNumber = {}
        geneCopyNumber['GCN0'] = []
        geneCopyNumber['GCN1'] = []
        geneCopyNumber['GCN2'] = []
        geneCopyNumber['GCN3'] = []
        geneCopyNumber['GCN4'] = []
        geneCopyNumber['GCN5+'] = []
        
        for marker in self.models:
            markerId = os.path.splitext(marker)[0]
            if marker in self.markerHits:
                if len(self.markerHits[marker]) >= 5:
                    geneCopyNumber['GCN5+'].append(markerId)
                else:
                    geneCopyNumber['GCN' + str(len(self.markerHits[marker]))].append(markerId)
            else:
                geneCopyNumber['GCN0'].append(markerId)
        
        return geneCopyNumber
          
    def getSummary(self, binMarkerSets, bIndividualMarkers, outputFormat=1):
        """Get dictionary containing information about bin."""
        summary = {}
        
        if outputFormat == 1:
            data, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            summary['marker lineage'] = bestMarkerSet.lineageStr
            summary['# markers'] = bestMarkerSet.numMarkers()
            summary['# marker sets'] = bestMarkerSet.numSets()
            summary['0'] = data[0]
            summary['1'] = data[1]
            summary['2'] = data[2]
            summary['3'] = data[3]
            summary['4'] = data[4]
            summary['5+'] = data[5]
            summary['Completeness'] = data[6]
            summary['Contamination'] = data[7]
        elif outputFormat == 2:
            data, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            summary['marker lineage'] = bestMarkerSet.lineageStr
            summary['# markers'] = bestMarkerSet.numMarkers()
            summary['# marker sets'] = bestMarkerSet.numSets()
            summary['0'] = data[0]
            summary['1'] = data[1]
            summary['2'] = data[2]
            summary['3'] = data[3]
            summary['4'] = data[4]
            summary['5+'] = data[5]
            summary['Completeness'] = data[6]
            summary['Contamination'] = data[7]
            summary.update(self.binStats)   
        elif outputFormat == 4:
            _, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            data = self.hitsToMarkerGene(bestMarkerSet)
            for marker, count in data.iteritems():
                summary[marker] = count

        elif outputFormat == 5:
            _, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            data = self.hitsToMarkerGene(bestMarkerSet)
            for marker, count in data.iteritems():
                summary[marker] = count

        elif outputFormat == 6:
            # tabular of bin_id, marker, contig_id
            for marker, hit_list in self.markerHits.items():
                summary[marker] = []
                for hit in hit_list:
                    summary[marker].append(hit.target_name)
                    
        elif outputFormat == 7:
            for marker, hit_list in self.markerHits.items():
                if len(hit_list) >= 2:
                    summary[marker] = []
                    for hit in hit_list:
                        summary[marker].append(hit.target_name)

        elif outputFormat == 8:
            # tabular - print only contigs that have more than one copy 
            # of the same marker on them
            contigs = defaultdict(dict)
            for marker, hit_list in self.markerHits.items():
                for hit in hit_list:
                    try:
                        contigs[hit.target_name][marker] += 1
                    except KeyError:
                        contigs[hit.target_name][marker] = 1
            
            for contig_name, marker_counts in contigs.items():
                for marker_name, marker_count in marker_counts.items():
                    if marker_count > 1:
                        if contig_name not in summary:
                            summary[contig_name] = {}
                            
                        summary[contig_name][marker_name] = marker_count
    
        elif outputFormat == 9:
            # tabular - print only position of marker genes
            genesWithMarkers = {}
            for marker, hit_list in self.markerHits.items():
                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(hit.target_name, []) + [hit]
                    
            for geneId, hits in genesWithMarkers.iteritems():
                summary[geneId] = {}
                for hit in hits:
                    summary[geneId][hit.query_accession] = summary[geneId].get(hit.query_accession, []) + [[hit.ali_from, hit.ali_to]]
                    
        elif outputFormat == 10:
            pass
        else:
            print("Unknown output format: ", outputFormat)
            
        return summary

    def printSummary(self, outputFormat, aai, binMarkerSets, bIndividualMarkers, coverageBinProfiles = None, table = None):
        """Print out information about bin."""
        if outputFormat == 1:
            data, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            row = "%s\t%s\t%d\t%d\t%s\t%0.2f\t%0.2f\t%0.2f" % (self.binId, bestMarkerSet.lineageStr,
                                                bestMarkerSet.numMarkers(), bestMarkerSet.numSets(),
                                                "\t".join([str(data[i]) for i in range(6)]),
                                                data[6],
                                                data[7],
                                                aai.aaiMeanBinHetero.get(self.binId, 0.0)
                                                )
            if table == None:
                print(row)
            else:  
                table.add_row([self.binId, bestMarkerSet.lineageStr, bestMarkerSet.numMarkers(), bestMarkerSet.numSets()] + data + [aai.aaiMeanBinHetero.get(self.binId, 0.0)])
        elif outputFormat == 2:
            data, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            
            if table == None:
                row = self.binId
                row += '\t%s\t%d\t%d' % (bestMarkerSet.lineageStr, bestMarkerSet.numMarkers(), bestMarkerSet.numSets())
                row += '\t%0.2f\t%0.2f\t%0.2f' % (data[6], data[7], aai.aaiMeanBinHetero.get(self.binId, 0.0))
                row += '\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (self.binStats['Genome size'], self.binStats['# scaffolds'], 
                                                 self.binStats['# contigs'], self.binStats['N50 (scaffolds)'], self.binStats['N50 (contigs)'], 
                                                 self.binStats['Longest scaffold'], self.binStats['Longest contig'])
                row += '\t%.1f\t%.2f' % (self.binStats['GC']*100, self.binStats['GC std']*100)
                row += '\t%.2f\t%d' % (self.binStats['Coding density'], self.binStats['# predicted ORFs'])
                row += '\t' + '\t'.join([str(data[i]) for i in xrange(6)])
                
                if coverageBinProfiles:
                    for _, coverageStats in coverageBinProfiles[self.binId].iteritems():
                        row += '\t%.2f\t%.2f' % (coverageStats[0], coverageStats[1])
            
                print(row)
            else:  
                row = [self.binId, bestMarkerSet.lineageStr, bestMarkerSet.numMarkers(), bestMarkerSet.numSets()]
                row.extend([data[6], data[7], aai.aaiMeanBinHetero.get(self.binId, 0.0)])
                row.extend([self.binStats['Genome size'], self.binStats['# scaffolds'], 
                                                 self.binStats['# contigs'], self.binStats['N50 (scaffolds)'], self.binStats['N50 (contigs)'], 
                                                 self.binStats['Longest scaffold'], self.binStats['Longest contig']])
                row.extend([self.binStats['GC']*100, self.binStats['GC std']*100])
                row.extend([self.binStats['Coding density'], self.binStats['# predicted ORFs']])
                row.extend(data[0:6])
                
                if coverageBinProfiles:
                    for _, coverageStats in coverageBinProfiles[self.binId].iteritems():
                        row.extend(coverageStats)
            
                table.add_row(row)
        elif outputFormat == 3:
            for ms in binMarkerSets.markerSetIter():
                data = self.geneCounts(ms, self.markerHits, bIndividualMarkers)
                row = "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%0.2f\t%0.2f\t%0.2f" % (self.binId, ms.UID, ms.lineageStr, ms.numGenomes,
                                                    ms.numMarkers(), ms.numSets(),
                                                    "\t".join([str(data[i]) for i in range(6)]),
                                                    data[6],
                                                    data[7],
                                                    aai.aaiMeanBinHetero.get(self.binId, 0.0)
                                                    )
                if table == None:
                    print(row)
                else:  
                    table.add_row([self.binId, ms.UID, ms.lineageStr, ms.numGenomes, ms.numMarkers(), ms.numSets()] + data + [aai.aaiMeanBinHetero.get(self.binId, 0.0)])
                    
        elif outputFormat == 4:
            _, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            data = self.hitsToMarkerGene(bestMarkerSet)
            print("--------------------")
            print(self.binId)
            
            markerHits = 0
            for marker, count in data.iteritems():
                print("%s\t%d" % (marker, count))
                if count >= 1:
                    markerHits += 1

            print("Total: %d / %d (%0.2f" % (markerHits,
                                              bestMarkerSet.numMarkers(),
                                              100*float(markerHits)/float(bestMarkerSet.numMarkers())
                                              )+"%)")
        elif outputFormat == 5:
            # matrix of bin vs marker counts
            _, bestMarkerSet = self.geneCountsForBestMarkerSet(binMarkerSets, bIndividualMarkers)
            data = self.hitsToMarkerGene(bestMarkerSet)
            
            columns = self.models.keys()
            
            rowStr = self.binId
            for marker in columns:
                count = 0
                try:
                    count = data[marker]
                except KeyError:
                    pass
                else:
                    rowStr += '\t' + str(count)
            print(rowStr)

        elif outputFormat == 6:
            # tabular of bin_id, marker, contig_id
            for marker, hit_list in self.markerHits.items():
                for hit in hit_list:
                    print(self.binId, marker, hit.target_name, sep='\t', end='\n')
                    
        elif outputFormat == 7:
            for marker, hitList in self.markerHits.items():
                if len(hitList) >= 2:
                    print(self.binId, marker, sep='\t', end='\t')
                    
                    scaffoldIds = []
                    for hit in hitList:
                        scaffoldIds.append(hit.target_name)
                        
                    print(','.join(sorted(scaffoldIds)), end='\n')

        elif outputFormat == 8:
            for marker, hitList in self.markerHits.items():
                if len(hitList) >= 2:
                    scaffoldsWithMultipleHits = set()
                    for i in xrange(0, len(hitList)):
                        scaffoldId = hitList[i].target_name[0:hitList[i].target_name.rfind('_')]
                        for j in xrange(i+1, len(hitList)):
                            if scaffoldId == hitList[j].target_name[0:hitList[j].target_name.rfind('_')]:
                                scaffoldsWithMultipleHits.add(hitList[i].target_name)
                                scaffoldsWithMultipleHits.add(hitList[j].target_name)

                    if len(scaffoldsWithMultipleHits) >= 2:
                        print(self.binId, marker, sep='\t', end='\t')
                        print(','.join(sorted(list(scaffoldsWithMultipleHits))), end='\n')
                    
        elif outputFormat == 9:
            # tabular - print only position of marker genes
            genesWithMarkers = {}
            for marker, hit_list in self.markerHits.items():
                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(hit.target_name, []) + [hit]
                    
            for geneId, hits in genesWithMarkers.iteritems():
                rowStr = self.binId + '\t' + geneId
                for hit in hits:
                    rowStr += '\t' + hit.ssion + ',' + str(hit.ali_from) + ',' + str(hit.ali_to)
                print(rowStr)
                    
        elif outputFormat == 10:
            markersInScaffold = {}
            for marker, hit_list in self.markerHits.items():
                for hit in hit_list:
                    scaffoldId = hit.target_name[0:hit.target_name.rfind('_')]
                    markersInScaffold[scaffoldId] = markersInScaffold.get(scaffoldId, []) + [marker]
            
            for scaffoldId, data in self.scaffoldStats.iteritems():
                print(scaffoldId, self.binId, str(data['Length']), str(data['# contigs']), 
                      '%.3f' % (data['GC']), str(data.get('# ORFs', 0)), 
                      '%.3f' % (float(data.get('Coding bases', 0)) / data['Total contig length']), 
                      sep='\t', end='\t')
                
                if scaffoldId in markersInScaffold:
                    markerStr = ','.join(sorted(markersInScaffold[scaffoldId]))
                    print(markerStr, end='\n')
                else:
                    print()
        else:
            self.logger.error("Unknown output format: %d", outputFormat)
