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
import uuid
import shutil
import ast
from collections import defaultdict
import logging

import defaultValues 
from common import reassignStdOut, restoreStdOut, checkFileExists

from hmmer import HMMERRunner, HMMERParser, makeOutputFNs
from hmmerModelParser import HmmModelParser

import prettytable

class ResultsParser():
    """Parse output of Prodigal+HMMER run and derived statistics."""
    def __init__(self):
        self.logger = logging.getLogger()
        
        # make the output file names
        self.results = {}
        self.models = {}
    
    def analyseResults(self,
                       directory,
                       hmmFile,
                       bIgnoreThresholds = False,
                       evalueThreshold = defaultValues.__CHECKM_DEFAULT_E_VAL__,
                       lengthThreshold = defaultValues.__CHECKM_DEFAULT_LENGTH__,
                       bSkipOrfCorrection = False,
                       ):
        """Parse the results in the output directory"""

        # parse the hmm file itself so we can determine the length of the queries
        modelParser = HmmModelParser(hmmFile)
        for model in modelParser.parse():
            self.models[model.name] = model
        
        # read bin and sequence stats into dictionaries
        binStats = self.parseBinStats(directory)
        seqStats = self.parseSeqStats(directory)

        # we expect directory to contain a collection of folders names after the original bins
        for folder in os.listdir(directory): 
            binFolder = os.path.join(directory, folder)
            if os.path.isdir(binFolder):
                # check if directory is a bin
                hmmerFile = os.path.join(binFolder, defaultValues.__CHECKM_DEFAULT_HMMER_TXT_OUT__)
                if not os.path.exists(hmmerFile):
                    continue
                
                resultsManager = ResultsManager(folder, bIgnoreThresholds, evalueThreshold, lengthThreshold, self.models, binStats[folder], seqStats[folder])
                       
                self.parseHmmerResults(hmmerFile, resultsManager, bSkipOrfCorrection)
                self.results[folder] = resultsManager
                
        # cache critical results to file
        self.writeBinStatsExt(directory)
        self.writeMarkerGeneStats(directory)

    def writeBinStatsExt(self, directory):
        binStatsExt = {}
        for binId in self.results:
            binStatsExt[binId] = self.results[binId].getSummary(outputFormat=1)
           
        binStatsExtFile = os.path.join(directory, 'storage', defaultValues.__CHECKM_DEFAULT_BIN_STATS_EXT_FILE__) 
        fout = open(binStatsExtFile, 'w')
        fout.write(str(binStatsExt))
        fout.close
        
    def writeMarkerGeneStats(self, directory):
        markerGenes = {}
        for binId in self.results:
            markerGenes[binId] = self.results[binId].getSummary(outputFormat=8)
           
        markerGenesFile = os.path.join(directory, 'storage', defaultValues.__CHECKM_DEFAULT_MARKER_GENE_STATS__)
        fout = open(markerGenesFile, 'w')
        fout.write(str(markerGenes))
        fout.close
                         
    def parseBinStats(self, resultsFolder):
        """Read bin statistics from file."""
        binStatsFile = os.path.join(resultsFolder, 'storage', defaultValues.__CHECKM_DEFAULT_BIN_STATS_FILE__)
        
        checkFileExists(binStatsFile)
        
        with open(binStatsFile, 'r') as f:
            s = f.read()
            binStats = ast.literal_eval(s)
            
        return binStats
    
    def parseBinStatsExt(self, resultsFolder):
        """Read bin statistics from file."""
        binStatsExtFile = os.path.join(resultsFolder, 'storage', defaultValues.__CHECKM_DEFAULT_BIN_STATS_EXT_FILE__)
        
        checkFileExists(binStatsExtFile)
        
        with open(binStatsExtFile, 'r') as f:
            s = f.read()
            binStatsExt = ast.literal_eval(s)
            
        return binStatsExt
    
    def parseMarkerGeneStats(self, resultsFolder):
        """Read bin statistics from file."""
        markerGeneStatsFile = os.path.join(resultsFolder, 'storage', defaultValues.__CHECKM_DEFAULT_MARKER_GENE_STATS__)
        
        checkFileExists(markerGeneStatsFile)
        
        with open(markerGeneStatsFile, 'r') as f:
            s = f.read()
            markerGeneStats = ast.literal_eval(s)
            
        return markerGeneStats
            
    def parseSeqStats(self, resultsFolder):
        """Read sequence statistics from file."""
        seqStatsFile = os.path.join(resultsFolder, 'storage', defaultValues.__CHECKM_DEFAULT_SEQ_STATS_FILE__)
        
        checkFileExists(seqStatsFile)
        
        with open(seqStatsFile, 'r') as f:
            s = f.read()
            seqStats = ast.literal_eval(s)
            
        return seqStats
            
    def parseHmmerResults(self, fileName, storage, bSkipOrfCorrection):
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
                    storage.addHit(hit)
                  
            if not bSkipOrfCorrection:
                storage.identifyOrfErrors()
                      
        except IOError as detail:
            sys.stderr.write(str(detail)+"\n")

    def getHeader(self, outputFormat):
        """Get header for requested output table."""
                
        # keep count of single, double, triple genes etc...
        if outputFormat == 1:
            header = ['Bin Id','0','1','2','3','4','5+','Completeness','Contamination']
        elif outputFormat == 2:
            header = ['Bin Id']
            header += ['Completeness','Contamination']
            header += ['Genome size (bp)', '# scaffolds', '# contigs', 'N50 (scaffolds)', 'N50 (contigs)', 'Longest scaffold (bp)', 'Longest contig (bp)']
            header += ['GC', 'GC std (scaffolds > 1Kbps)']
            header += ['Coding density (translation table 11)', '# predicted ORFs']
            header += ['0','1','2','3','4','5+']

            if defaultValues.__CHECKM_DEFAULT_MIN_SEQ_LEN_GC_STD__ != 1000:
                self.logger.error('  [Error] Labeling error: GC std (scaffolds > 1Kbps)')
                sys.exit() 
        elif outputFormat == 3:
            header = []
        elif outputFormat == 4:
            header = [''] + self.models.keys()
        elif outputFormat == 5:
            header = ['Bin Id','Marker Id','Scaffold Id']
        elif outputFormat == 6:
            header = ['Bin Id','Marker Id','Scaffold Ids']
        elif outputFormat == 7:
            header = ['Bin Id','Scaffold Id','{Marker Id, Count}']
        elif outputFormat == 8:
            header = ['Bin Id','Gene Id','{Marker Id, Start position, End position}']
        elif outputFormat == 9:
            header = ['Scaffold Id','Bin Id','Length','# contigs','GC','# ORFs','Coding density','Marker Ids']
            
        return header
        
    def printSummary(self, outputFormat, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)
        
        print('')
        
        prettyTableFormats = [1, 2]      
          
        header = self.getHeader(outputFormat) 
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
            self.results[binId].printSummary(outputFormat, pTable)
            
        if not bTabTable :  
            print(pTable.get_string())
        
        print('') 
            
        # restore stdout   
        restoreStdOut(outFile, oldStdOut)     

class ResultsManager():
    """Store all the results for a single bin"""
    def __init__(self, binId, bIgnoreThresholds, evalueThreshold, lengthThreshold, models=None, binStats=None, scaffoldStats=None):
        self.binId = binId
        self.markers = {}
        self.bIgnoreThresholds = bIgnoreThresholds
        self.evalueThreshold = evalueThreshold
        self.lengthThreshold = lengthThreshold
        self.models = models
        self.binStats = binStats
        self.scaffoldStats = scaffoldStats
    
    def vetHit(self, hit):
        """Check if hit meets required thresholds."""
        
        model = self.models[hit.query_name]
        
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
            try:
                # retain only the best domain hit for a given marker to a specific ORF
                previousHitToORF = None
                for h in self.markers[hit.query_name]:
                    if h.target_name == hit.target_name:
                        previousHitToORF = h
                        break
                    
                if not previousHitToORF:
                    self.markers[hit.query_name].append(hit)
                else:
                    if previousHitToORF.dom_score < hit.dom_score:
                        self.markers[hit.query_name].append(hit)
                        self.markers[hit.query_name].remove(previousHitToORF)
                    
            except KeyError:
                self.markers[hit.query_name] = [hit]
                
    def identifyOrfErrors(self):
        """Identify ORF errors affecting marker genes."""
 
        # check for adjacent ORFs with hits to the same marker gene
        for markerId, hits in self.markers.iteritems():
            
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
                        newHit.target_name = '|'.join(sorted([orfI, orfJ]))
                        
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
                    
            self.markers[markerId] = hits
                
    def calculateMarkers(self, verbose=False):
        """Returns an object containing summary information 
           When verbose is False a list is returned containing the counts of
           markers for the bin as well as the total completeness and
           contamination.  If the verbose flag is set to true, a dict is
           returned containing each marker as the key and the value as the
           count.
        """
        if verbose:
            ret = dict()
            for marker in self.models:
                try:
                    ret[marker] = len(self.markers[marker])
                except KeyError:
                    ret[marker] = 0
            return ret
        else:
            geneCounts = [0]*6
            multiCopyCount = 0
            for marker in self.models:
                # we need to limit it form 0 to 5+
                try:
                    if len(self.markers[marker]) > 5:
                        markerCount = 5
                    else:
                        markerCount = len(self.markers[marker])
                        
                    multiCopyCount += (len(self.markers[marker]) - 1)
                except KeyError:
                    markerCount = 0
                
                geneCounts[markerCount] += 1
                
            percComp = 100 * float(len(self.markers)) / float(len(self.models))
            percCont = 100 * float(multiCopyCount) / float(len(self.models))  
            
            geneCounts.append(percComp)
            geneCounts.append(percCont)
            
            return geneCounts
        
    def getSummary(self, outputFormat=1):
        """Get dictionary containing information about bin."""
        summary = {}
        
        if outputFormat == 1:
            data = self.calculateMarkers(verbose=False)
            summary['0'] = data[0]
            summary['1'] = data[1]
            summary['2'] = data[2]
            summary['3'] = data[3]
            summary['4'] = data[4]
            summary['5+'] = data[5]
            summary['Completeness'] = data[6]
            summary['Contamination'] = data[7]
        elif outputFormat == 1:
            data = self.calculateMarkers(verbose=False)
            summary['0'] = data[0]
            summary['1'] = data[1]
            summary['2'] = data[2]
            summary['3'] = data[3]
            summary['4'] = data[4]
            summary['5+'] = data[5]
            summary['Completeness'] = data[6]
            summary['Contamination'] = data[7]
            summary.update(self.binStats)   
        elif outputFormat == 3:
            data = self.calculateMarkers(verbose=True)
            for marker,count in data.iteritems():
                summary[marker] = count

        elif outputFormat == 4:
            data = self.calculateMarkers(verbose=True)
            for marker,count in data.iteritems():
                summary[marker] = count

        elif outputFormat == 5:
            # tabular of bin_id, marker, contig_id
            for marker, hit_list in self.markers.items():
                summary[marker] = []
                for hit in hit_list:
                    summary[marker].append(hit.target_name)
                    
        elif outputFormat == 6:
            for marker, hit_list in self.markers.items():
                if len(hit_list) >= 2:
                    summary[marker] = []
                    for hit in hit_list:
                        summary[marker].append(hit.target_name)

        elif outputFormat == 7:
            # tabular - print only contigs that have more than one copy 
            # of the same marker on them
            contigs = defaultdict(dict)
            for marker, hit_list in self.markers.items():
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
    
        elif outputFormat == 8:
            # tabular - print only position of marker genes
            genesWithMarkers = {}
            for marker, hit_list in self.markers.items():
                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(hit.target_name, []) + [hit]
                    
            for geneId, hits in genesWithMarkers.iteritems():
                summary[geneId] = {}
                for hit in hits:
                    summary[geneId][hit.query_name] = summary[geneId].get(hit.query_name, []) + [[hit.ali_from, hit.ali_to]]
                    
        elif outputFormat == 9:
            pass
        else:
            print("Unknown output format: ", outputFormat)
            
        return summary

    def printSummary(self, outputFormat, table = None):
        """Print out information about bin."""
        if outputFormat == 1:
            data = self.calculateMarkers(verbose=False)
            row = "%s\t%s\t%0.2f\t%0.2f" % (self.binId,
                                                "\t".join([str(data[i]) for i in range(6)]),
                                                data[6],
                                                data[7]
                                                )
            if table == None:
                print(row)
            else:  
                table.add_row([self.binId] + data)
        elif outputFormat == 2:
            data = self.calculateMarkers(verbose=False)
            
            if table == None:
                row = self.binId
                row += '\t%0.2f\t%0.2f' % (data[6], data[7])
                row += '\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (self.binStats['Genome size'], self.binStats['# scaffolds'], 
                                                 self.binStats['# contigs'], self.binStats['N50 (scaffolds)'], self.binStats['N50 (contigs)'], 
                                                 self.binStats['Longest scaffold'], self.binStats['Longest contig'])
                row += '\t%.1f\t%.2f' % (self.binStats['GC']*100, self.binStats['GC std']*100)
                row += '\t%.2f\t%d' % (self.binStats['Coding density'], self.binStats['# predicted ORFs'])
                row += '\t' + '\t'.join([str(data[i]) for i in xrange(6)])
            
                print(row)
            else:  
                row = [self.binId]
                row.extend([data[6], data[7]])
                row.extend([self.binStats['Genome size'], self.binStats['# scaffolds'], 
                                                 self.binStats['# contigs'], self.binStats['N50 (scaffolds)'], self.binStats['N50 (contigs)'], 
                                                 self.binStats['Longest scaffold'], self.binStats['Longest contig']])
                row.extend([self.binStats['GC']*100, self.binStats['GC std']*100])
                row.extend([self.binStats['Coding density'], self.binStats['# predicted ORFs']])
                row.extend(data[0:6])
            
                table.add_row(row)
                    
        elif outputFormat == 3:
            data = self.calculateMarkers(verbose=True)
            print("--------------------")
            print(self.binId)
            for marker,count in data.iteritems():
                print("%s\t%d" % (marker, count))

            print("TOTAL:\t%d / %d (%0.2f" % (len(self.markers),
                                              len(self.models),
                                              100*float(len(self.markers))/float(len(self.models))
                                              )+"%)")
        elif outputFormat == 4:
            # matrix of bin vs marker counts
            data = self.calculateMarkers(verbose=True)
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

        elif outputFormat == 5:
            # tabular of bin_id, marker, contig_id
            for marker, hit_list in self.markers.items():
                for hit in hit_list:
                    print(self.binId, marker, hit.target_name, sep='\t', end='\n')
                    
        elif outputFormat == 6:
            for marker, hitList in self.markers.items():
                if len(hitList) >= 2:
                    print(self.binId, marker, sep='\t', end='\t')
                    
                    scaffoldIds = []
                    for hit in hitList:
                        scaffoldIds.append(hit.target_name)
                        
                    print(','.join(sorted(scaffoldIds)), end='\n')

        elif outputFormat == 7:
            for marker, hitList in self.markers.items():
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
                    
        elif outputFormat == 8:
            # tabular - print only position of marker genes
            genesWithMarkers = {}
            for marker, hit_list in self.markers.items():
                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(hit.target_name, []) + [hit]
                    
            for geneId, hits in genesWithMarkers.iteritems():
                rowStr = self.binId + '\t' + geneId
                for hit in hits:
                    rowStr += '\t' + hit.query_name + ',' + str(hit.ali_from) + ',' + str(hit.ali_to)
                print(rowStr)
                    
        elif outputFormat == 9:
            markersInScaffold = {}
            for marker, hit_list in self.markers.items():
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
            print("Unknown output format: ", outputFormat)

class HMMAligner:
    def __init__(self, individualFile=False, includeConsensus=True, outputFormat='PSIBLAST'):
        # make output file names
        self.individualFile = individualFile
        self.includeConsensus = includeConsensus
        self.outputFormat = outputFormat
        (self.txtOut, self.hmmOut) = makeOutputFNs()
        (_, self.hmmAlign) = makeOutputFNs(mode='align')
        
    def makeAlignments(self,
                       directory,
                       hmm,
                       evalueThreshold = defaultValues.__CHECKM_DEFAULT_E_VAL__,
                       lengthThreshold = defaultValues.__CHECKM_DEFAULT_LENGTH__,
                       bestHit=False,
                       generateModelFiles=True
                       ):
        """main wrapper for making hmm alignments
           first parse through the hmm txt and work out all the 
           hits that were made
           we expect directory to contain a collection of folders
           names after the original bins
           we need to use this guy to vet the hits
         """
        models = {}
        model_parser = HmmModelParser(hmm)
        for model in model_parser.parse():
            models[model.name] = model

        resultsManager = ResultsManager('', lengthThreshold, evalueThreshold, models)
        hit_lookup = {}
        unique_hits = {}
        for folder in os.listdir(directory):
            if not os.path.isdir(os.path.join(directory, folder)):
                continue
            
            # we can now build the hmmer_file_name
            hmmer_file_name = os.path.join(directory, folder, self.txtOut)
            hit_lookup[folder] = {}
            with open(hmmer_file_name, 'r') as hmmerHandle:
                try:
                    HP = HMMERParser(hmmerHandle)
                except:
                    print("Error opening HMM file: ", hmmer_file_name)
                    raise
                
                while True:
                    hit = HP.next()
                    if hit is None:
                        break
                    if resultsManager.vetHit(hit): # we only want good hits
                        try:
                            if bestHit:
                                if hit_lookup[folder][hit.query_name][0].full_e_value > hit.full_e_value:
                                    hit_lookup[folder][hit.query_name] = [hit]
                            else:
                                hit_lookup[folder][hit.query_name].append(hit)
                        except KeyError:
                            hit_lookup[folder][hit.query_name] = [hit]
                            
                        try:
                            unique_hits[hit.query_name] += 1
                        except KeyError:
                            unique_hits[hit.query_name] = 1
        
        # make temporary files and write the extracted hmms
        single_hmms = {}
        single_hmm_consensi = {}
        HA = HMMERRunner(mode='align')
        if generateModelFiles:
            HF = HMMERRunner(mode='fetch')
            self.makeAlignmentModels(HF, hmm, unique_hits.keys(), single_hmms, single_hmm_consensi)
        
        # now we need to go through each of the bins and extract the contigs which we need
        for folder in hit_lookup:
            genes_file_name = os.path.join(directory, folder,
                    defaultValues.__CHECKM_DEFAULT_PRODIGAL_AA__) # the seqs to align
            if 0 < len(hit_lookup[folder]):
                self.alignBin(HA, hit_lookup[folder], single_hmms,
                        single_hmm_consensi, genes_file_name, directory, folder)
            else:
                self.logger.info("Skipping alignment for", folder, ": no hits found")
                
        # remove the tmp files
        for file_name in single_hmms:
            os.remove(single_hmms[file_name])
        for file_name in single_hmm_consensi:
            os.remove(single_hmm_consensi[file_name])
            
    def alignBin(self, HA, hits, hmms, consensi, fasta,
            directory, folder):
        # hits is a hash of HMM names to hit objects
        """Make alignments for this bin"""
        final_alignment_file = os.path.join(directory, folder, self.hmmAlign) # where to write the whole shebang to
        # a hash of sequence names and sequences
        bin_genes = {}
        with open(fasta, 'r') as fh:
            for cid, seq, _ in self.readfq(fh):
                bin_genes[cid] = seq

        align_tmp_file = os.path.join('/tmp', str(uuid.uuid4()))
        for hit_name in hits:
            align_seq_file = open(os.path.join('/tmp', str(uuid.uuid4())), 'w')

            if self.includeConsensus:
                with open(consensi[hit_name]) as fh:
                    for line in fh:
                        align_seq_file.write(line)

            for hmm_hit in hits[hit_name]:
                start_location = 0
                end_location = len(bin_genes[hmm_hit.target_name])
                
                if hmm_hit.ali_to+100 < end_location:
                    end_location = hmm_hit.ali_to+100 
                
                if hmm_hit.ali_from-100 > start_location:
                    start_location = hmm_hit.ali_from-100 

                align_seq_file.write(">"+hmm_hit.target_name+"\n")
                align_seq_file.write(bin_genes[hmm_hit.target_name][start_location:end_location])
                align_seq_file.write("\n")
            
            align_seq_file.close()
            
            if self.individualFile:
                out_file = os.path.join(directory, folder,hit_name+"_out.align" )
                HA.align(hmms[hit_name], align_seq_file.name, out_file,
                        writeMode='>', outputFormat=self.outputFormat, allcol=False)
                
            else:
                HA.align(hmms[hit_name], align_seq_file.name, align_tmp_file,
                        writeMode='>>', outputFormat=self.outputFormat, allcol=False)
                with open(align_tmp_file, 'a') as fh:
                    fh.write("\n###########################################\n\n")
            
            os.remove(align_seq_file.name)

        if not self.individualFile:
            shutil.copy(align_tmp_file, final_alignment_file)
            os.remove(align_tmp_file)
        
    def makeAlignmentModels(self, HF, hmm, keyz, hmms, consensi):
        """Make a bunch of temporary hmm files which we'll use when aligning"""
        for key in keyz:
            fetch_filename = os.path.join('/tmp', str(uuid.uuid4()))
            emit_filename = os.path.join('/tmp', str(uuid.uuid4()))
            hmms[key] = fetch_filename
            consensi[key] = emit_filename
            HF.fetch(hmm, key, fetch_filename, emit_filename)
            
    def readfq(self, fp): # this is a generator function
        """https://github.com/lh3"""
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break