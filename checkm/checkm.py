###############################################################################
#
# checkm.py - wraps coarse workflows
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
import logging
from collections import defaultdict

import numpy as np

import defaultValues

from timeKeeper import TimeKeeper
from markerSets import MarkerSetParser
from resultsParser import ResultsParser
from hmmerAligner import HmmerAligner
from markerGeneFinder import MarkerGeneFinder
from pplacer import PplacerRunner
from treeParser import TreeParser
from taxonParser import TaxonParser
from aminoAcidIdentity import AminoAcidIdentity
from binComparer import BinComparer
from binStatistics import BinStatistics
from coverage import Coverage
from coverageWindows import CoverageWindows
from genomicSignatures import GenomicSignatures
from unbinned import Unbinned
from merger import Merger
from profile import Profile
from binTools import BinTools
from ssuFinder import SSU_Finder
from PCA import PCA
from testing import Testing
from common import makeSurePathExists, checkFileExists, binIdFromFilename, reassignStdOut, restoreStdOut, getBinIdsFromOutDir, checkDirExists

from plot.gcPlots import GcPlots
from plot.codingDensityPlots import CodingDensityPlots
from plot.tetraDistPlots import TetraDistPlots
from plot.distributionPlots import DistributionPlots
from plot.nxPlot import NxPlot
from plot.cumulativeLengthPlot import CumulativeLengthPlot
from plot.lengthHistogram import LengthHistogram
from plot.markerGenePosPlot import MarkerGenePosPlot
from plot.parallelCoordPlot import ParallelCoordPlot
from plot.binQAPlot import BinQAPlot
from plot.pcaPlot import PcaPlot
from plot.gcBiasPlots import GcBiasPlot

class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()   
    
    def binFiles(self, binFolder, binExtension):
        binFiles = []
        if binFolder is not None:
            all_files = os.listdir(binFolder)
            for f in all_files:
                if f.endswith(binExtension):
                    binFiles.append(os.path.join(binFolder, f))
                    
        if not binFiles:
            self.logger.error("  [Error] No bins found. Check the extension (-x) used to identify bins.")
            sys.exit()
        
        return sorted(binFiles)
    
    def tree(self, options):
        """Tree command"""   
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tree] Placing bins in genome tree.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
                
        binFiles = self.binFiles(options.bin_folder, options.extension)   
        
        # setup directory structure
        makeSurePathExists(options.out_folder)
        makeSurePathExists(os.path.join(options.out_folder, 'bins')) 
        makeSurePathExists(os.path.join(options.out_folder, 'storage')) 
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'hmms'))
        
        # find phylogenetically informative genes in genome bins 
        phyloHMMs = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'hmms', 'phylo.hmm')                              
        mgf = MarkerGeneFinder(options.threads)
        mgf.find(binFiles, options.out_folder, defaultValues.HMMER_TABLE_PHYLO_OUT, defaultValues.HMMER_PHYLO_OUT, phyloHMMs)
        
        # calculate statistics for each genome bin
        self.logger.info('')
        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.out_folder, defaultValues.BIN_STATS_PHYLO_OUT, defaultValues.SEQ_STATS_PHYLO_OUT)
        
        # set HMM file for each bin
        binIdToHmmModelFile = {}
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)
            binIdToHmmModelFile[binId] = phyloHMMs
        
        # align identified marker genes
        self.logger.info('')
        HA = HmmerAligner(options.threads)
        HA.makeAlignmentToCommonMarkers(options.out_folder,
                                          defaultValues.HMMER_TABLE_PHYLO_OUT,
                                          binIdToHmmModelFile,
                                          False,
                                          defaultValues.E_VAL,
                                          defaultValues.LENGTH,
                                          os.path.join(options.out_folder, 'storage', 'tree')
                                          )
        
        # place bins into genome tree
        self.logger.info('')
        pplacer = PplacerRunner(options.threads)
        pplacer.run(binFiles, options.out_folder)
                
        self.timeKeeper.printTimeStamp()
        
    def treeQA(self, options):
        """QA command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tree_qa] Assessing phylogenetic markers found in each bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.tree_folder)
        
        # set HMM file for each bin
        phyloHMMs = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'hmms', 'phylo.hmm') 
        markerSetParser = MarkerSetParser()
        binIdToHmmModelFile = markerSetParser.createHmmModelFiles(options.tree_folder, getBinIdsFromOutDir(options.tree_folder), phyloHMMs)
        
        # calculate marker gene statistics
        RP = ResultsParser()
        RP.analyseResults(options.tree_folder, 
                          defaultValues.BIN_STATS_PHYLO_OUT, 
                          defaultValues.SEQ_STATS_PHYLO_OUT, 
                          defaultValues.HMMER_TABLE_PHYLO_OUT, 
                          binIdToHmmModelFile)

        # determine taxonomy of each bin
        treeParser = TreeParser()
        treeParser.printSummary(options.out_format, options.tree_folder, RP, options.bTabTable, options.file)
         
        if options.file != '':
            self.logger.info('  QA information written to: ' + options.file)
            
        self.timeKeeper.printTimeStamp()

    def lineageSet(self, options, db=None):
        """Lineage set command"""   
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - lineage_set] Inferring lineage-specific marker sets.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
                
        checkDirExists(options.tree_folder) 

        treeParser = TreeParser()
        treeParser.getBinMarkerSets(options.tree_folder, options.marker_file, 
                                    options.num_genomes_markers, options.num_genomes_refine, 
                                    options.bootstrap, options. bNoLineageSpecificRefinement,
                                    options.bRequireTaxonomy)
        
        self.logger.info('')
        self.logger.info('  Marker set written to: ' + options.marker_file)
        
        self.timeKeeper.printTimeStamp()
        
    def taxonList(self, options, db=None):
        """Lineage set command"""   
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - taxon_list] Listing available taxonomic-specific marker sets.')
        self.logger.info('*******************************************************************************')
                
        taxonParser = TaxonParser()
        taxonParser.list(options.rank)
        
        self.timeKeeper.printTimeStamp()
        
    def taxonSet(self, options, db=None):
        """Taxon set command"""   
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - taxon_set] Inferring taxonomic-specific marker set.')
        self.logger.info('*******************************************************************************')          
        self.logger.info('')
        
        taxonParser = TaxonParser()
        bValidSet = taxonParser.markerSet(options.rank, options.taxon, options.marker_file)
        
        if bValidSet:
            self.logger.info('')
            self.logger.info('  Marker set written to: ' + options.marker_file)
        
        self.timeKeeper.printTimeStamp()
        
    def analyze(self, options, db=None):
        """Analyze command"""   
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - analyze] Identifying marker genes in bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
                
        binFiles = self.binFiles(options.bin_folder, options.extension)   

        # setup directory structure
        makeSurePathExists(options.out_folder)
        makeSurePathExists(os.path.join(options.out_folder, 'bins')) 
        makeSurePathExists(os.path.join(options.out_folder, 'storage')) 
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'aai_qa'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'hmms'))
        
        # find marker genes in genome bins    
        mgf = MarkerGeneFinder(options.threads)
        binIdToHmmModelFile = mgf.find(binFiles, options.out_folder, defaultValues.HMMER_TABLE_OUT, defaultValues.HMMER_OUT, options.marker_file)
        
        markerSetParser = MarkerSetParser(options.threads)
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.out_folder, getBinIdsFromOutDir(options.out_folder), options.marker_file)
        
        self.timeKeeper.printTimeStamp()
        
        # align marker genes with multiple hits within a bin
        self.logger.info('')
        HA = HmmerAligner(options.threads)
        HA.makeAlignmentsOfMultipleHits(options.out_folder,
                                          defaultValues.HMMER_TABLE_OUT,
                                          binIdToHmmModelFile,
                                          binIdToBinMarkerSets,
                                          False,
                                          defaultValues.E_VAL,
                                          defaultValues.LENGTH,
                                          os.path.join(options.out_folder, 'storage', 'aai_qa')
                                          )
        
        self.timeKeeper.printTimeStamp()
        
        # calculate statistics for each genome bin
        self.logger.info('')
        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.out_folder, defaultValues.BIN_STATS_OUT, defaultValues.SEQ_STATS_OUT)
        
        self.timeKeeper.printTimeStamp()

    def qa(self, options):
        """QA command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - qa] Tabulating genome statistics.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.analyze_folder)

        # calculate AAI between marks with multiple hits in a single bin
        aai = AminoAcidIdentity()
        aai.run(options.aai_strain, options.analyze_folder)
        
        # get HMM file for each bin
        markerSetParser = MarkerSetParser(options.threads)
        binIdToHmmModelFile = markerSetParser.createHmmModelFiles(options.analyze_folder, getBinIdsFromOutDir(options.analyze_folder), options.marker_file)
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.analyze_folder, getBinIdsFromOutDir(options.analyze_folder), options.marker_file)

        # get results for each bin
        RP = ResultsParser()
        RP.analyseResults(options.analyze_folder, 
                          defaultValues.BIN_STATS_OUT, 
                          defaultValues.SEQ_STATS_OUT,
                          defaultValues.HMMER_TABLE_OUT,
                          binIdToHmmModelFile,
                          bIgnoreThresholds = options.bIgnoreThresholds,
                          evalueThreshold = options.e_value,
                          lengthThreshold = options.length,
                          bSkipOrfCorrection = options.bSkipOrfCorrection
                          )
        
        self.logger.info('')
        RP.printSummary(options.out_format, aai, binIdToBinMarkerSets, options.bIndividualMarkers, options.coverage_file, options.bTabTable, options.file)
        RP.cacheResults(options.analyze_folder, binIdToBinMarkerSets, options.bIndividualMarkers)
                        
        if options.file != '':
            self.logger.info('  QA information written to: ' + options.file)
            
        self.timeKeeper.printTimeStamp()
                
    def gcPlot(self, options):
        """GC plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - gc_plot] Creating GC histogram and delta-GC plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
    
        plots = GcPlots(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
            
            plots.plot(f, options.distributions)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.gc_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def codingDensityPlot(self, options):
        """Coding density plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - coding_plot] Creating coding density histogram and delta-CD plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
            
        plots = CodingDensityPlots(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting coding density plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
            
            plots.plot(f, options.distributions)
            
            binId = binIdFromFilename(f) 
            outputFile = os.path.join(options.plot_folder, binId) + '.coding_density_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def tetraDistPlot(self, options):
        """Tetranucleotide distance plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra_plot] Creating tetra-distance histogram and delta-TD plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
        
        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)
            
        plots = TetraDistPlots(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting tetranuclotide distance plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
                
            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.tetra_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def distributionPlots(self, options):
        """Reference distribution plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - dist_plot] Creating GC, CD, and TD distribution plots.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
        
        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)
            
        plots = DistributionPlots(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting reference distribution plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
                
            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.ref_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def tetraPcaPlot(self, options):
        """PCA plot of tetranucleotide signatures"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra_pca] Creating PCA plot of tetranucleotide signatures.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
        
        self.logger.info('  Computing PCA of tetranuclotide signatures.\n')
        pca = PCA()
        seqIds, pc, variance = pca.pcaFile(options.tetra_profile, fraction=1.0, bCenter=True, bScale=False)
            
        plots = PcaPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting PCA of tetranuclotide signatures for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
                    
            plots.plot(f, seqIds, pc, variance)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.tetra_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def coveragePcaPlot(self, options):
        """PCA plot of coverage profiles"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - cov_pca] Creating PCA plot of coverage profiles.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        checkFileExists(options.coverage_file)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
        
        coverage = Coverage(threads=1)
        coverageStats = coverage.parseCoverage(options.coverage_file)
        
        seqIds = []
        coverageProfiles = []
        for binId, seqDict in coverageStats.iteritems():
            for seqId, bamDict in seqDict.iteritems():
                seqIds.append(seqId)
                
                coverages = []
                for _, coverage in bamDict.iteritems():  
                    coverages.append(coverage)

                coverageProfiles.append(coverages)
                
        coverageProfiles = np.array(coverageProfiles)
        if coverageProfiles.shape[1] < 2:
            self.logger.error('  [Error] Coverage profile is 1 dimensional. PCA requires at least 2 dimensions.')
            sys.exit()
        
        self.logger.info('  Computing PCA of coverage profiles.\n')
        pca = PCA()
        pc, variance = pca.pcaMatrix(coverageProfiles, fraction=1.0, bCenter=True, bScale=False)
            
        plots = PcaPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting PCA of coverage profiles for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
                    
            plots.plot(f, seqIds, pc, variance)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.cov_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def gcBiasPlot(self, options):
        """GC bias plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - gc_bias_plot] Plotting bin coverage as a function of GC.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension) 
        
        coverageWindows = CoverageWindows(options.threads)
        coverageProfile = coverageWindows.run(binFiles, options.bam_file, options.all_reads, options.min_align, options.max_edit_dist, options.window_size)
        
        plots = GcBiasPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            self.logger.info('  Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1
            
            plots.plot(f, coverageProfile)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.gc_bias_plot.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
            
    def nxPlot(self, options):
        """Nx-plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - nx_plot] Creating Nx-plots.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)  
            
        nx = NxPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting Nx-plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            nx.plot(f)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.nx_plot.' + options.image_type
            nx.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def cumulativeLengthPlot(self, options):
        """Cumulative sequence length plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - len_plot] Creating cumulative sequence length plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)
            
        plot = CumulativeLengthPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting cumulative sequence length plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            plot.plot(f)
             
            outputFile = os.path.join(options.plot_folder, binId) + '.seq_len_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def lengthHistogram(self, options):
        """Sequence length histogram command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - len_hist] Creating sequence length histogram.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
            
        binFiles = self.binFiles(options.bin_folder, options.extension)
            
        plot = LengthHistogram(options)
        filesProcessed = 1
        for f in binFiles:  
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting sequence length histogram for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            plot.plot(f)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.len_hist.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def markerPlot(self, options):
        """Marker gene position plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - marker_plot] Creating marker gene position plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
                        
        # generate plot for each bin    
        binFiles = self.binFiles(options.bin_folder, options.extension)
        
        resultsParser = ResultsParser()
        markerGeneStats = resultsParser.parseMarkerGeneStats(options.out_folder)
        binStats = resultsParser.parseBinStatsExt(options.out_folder)
            
        plot = MarkerGenePosPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting marker gene position plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
              
            bPlotted = plot.plot(f, markerGeneStats[binId], binStats[binId])
            
            if bPlotted:
                outputFile = os.path.join(options.plot_folder, binId) + '.marker_pos_plot.' + options.image_type
                plot.savePlot(outputFile, dpi=options.dpi)
                self.logger.info('    Plot written to: ' + outputFile)
            else:
                self.logger.info('    No marker genes found in bin.')
            
        self.timeKeeper.printTimeStamp()

    def parallelCoordPlot(self, options):
        """Parallel coordinate plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - par_plot] Creating parallel coordinate plot of GC and coverage.')
        self.logger.info('*******************************************************************************')  
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
        checkFileExists(options.coverage_file)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        # read sequence stats file
        resultsParser = ResultsParser()
        seqStats = resultsParser.parseSeqStats(options.out_folder, defaultValues.SEQ_STATS_OUT)

        # read coverage stats file
        coverage = Coverage(threads = 1)
        coverageStats = coverage.parseCoverage(options.coverage_file)
        
        # create plot for each bin
        plot = ParallelCoordPlot(options)
        filesProcessed = 1
        for f in binFiles:  
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting marker gene position plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
                
            plot.plot(binId, seqStats, coverageStats)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.paralel_coord_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def binQAPlot(self, options):
        """Bin QA plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - bin_qa_plot] Creating bar plot of bin quality.')
        self.logger.info('*******************************************************************************')  
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)
        
        # read sequence stats file
        resultsParser = ResultsParser()
        binStatsExt = resultsParser.parseBinStatsExt(options.out_folder)
  
        # create plot for each bin
        plot = BinQAPlot(options)
        if not options.bIgnoreHetero:
            aai = AminoAcidIdentity()
            aai.run(options.aai_strain, options.out_folder)
            plot.plot(binFiles, binStatsExt, options.bIgnoreHetero, aai.aaiHetero)
        else:
            plot.plot(binFiles, binStatsExt, options.bIgnoreHetero, None)
        
        outputFile = os.path.join(options.plot_folder, 'bin_qa_plot.' + options.image_type)
        plot.savePlot(outputFile, dpi=options.dpi)
        self.logger.info('')
        self.logger.info('  Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def unbinned(self, options):
        """Unbinned Command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - unbinned] Identify unbinned sequences.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
  
        binFiles = self.binFiles(options.bin_folder, options.extension)
  
        unbinned = Unbinned()
        unbinned.run(binFiles, options.seq_file, options.output_seq_file, options.output_stats_file, options.min_seq_len)
        
        self.logger.info('')
        self.logger.info('  Unbinned sequences written to: ' + options.output_seq_file)
        self.logger.info('  Unbinned sequences statistics written to: ' + options.output_stats_file)
            
        self.timeKeeper.printTimeStamp()
        
    def coverage(self, options):
        """Coverage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - coverage] Calculating coverage of sequences.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
    
        binFiles = self.binFiles(options.bin_folder, options.extension) 
        
        coverage = Coverage(options.threads)
        coverage.run(binFiles, options.bam_files, options.output_file, options.all_reads, options.min_align, options.max_edit_dist)
        
        self.logger.info('  Coverage information written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def tetraSignatures(self, options):
        """Tetranucleotide signature command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra] Calculating tetranucleotide signature of sequences.')
        self.logger.info('*******************************************************************************')
  
        checkFileExists(options.seq_file)
        
        self.logger.info('')
        tetraSig = GenomicSignatures(4, options.threads)
        tetraSig.calculate(options.seq_file, options.output_file)
        
        self.logger.info('  Tetranucletoide signatures written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def profile(self, options):
        """Profile command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - profile] Calculating percentage of reads mapped to each bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkFileExists(options.coverage_file)
  
        profile = Profile()
        profile.run(options.coverage_file, options.file, options.bTabTable)
        
        if options.file != '':
            self.logger.info('  Profile information written to: ' + options.file)
          
        self.timeKeeper.printTimeStamp()
        
    def merge(self, options):
        """Merge command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - merge] Identifying bins with complementary sets of marker genes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        
        # setup directory structure
        makeSurePathExists(options.out_folder)
        makeSurePathExists(os.path.join(options.out_folder, 'bins')) 
        makeSurePathExists(os.path.join(options.out_folder, 'storage')) 
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'aai_qa'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'hmms'))
        
        binFiles = self.binFiles(options.bin_folder, options.extension)
        
        binIds = []
        for binFile in binFiles:
            binIds.append(binIdFromFilename(binFile))
        
        # find marker genes in genome bins                         
        mgf = MarkerGeneFinder(options.threads)
        binIdToHmmModelFile = mgf.find(binFiles, options.out_folder, "merger.table.txt", "merger.hmmer3", options.marker_file)
        
        # get HMM file for each bin
        markerSetParser = MarkerSetParser()
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.out_folder, binIds, options.marker_file)

        # compare markers found in each bin
        self.logger.info('')
        merger = Merger()
        outputFile = merger.run(binFiles, options.out_folder, "merger.table.txt", binIdToHmmModelFile, binIdToBinMarkerSets,
                                options.delta_comp, options.delta_cont, options.merged_comp, options.merged_cont)
        
        self.logger.info('\n  Merger information written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
            
    def outliers(self, options):
        """Outlier command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - outlier] Identifying outliers in bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        checkDirExists(options.bin_folder)
        checkFileExists(options.tetra_profile)
        
        binFiles = self.binFiles(options.bin_folder, options.extension)
        
        binTools = BinTools()
        binTools.identifyOutliers(options.out_folder, binFiles, options.tetra_profile, options.distribution, options.report_type, options.output_file)

        self.logger.info('\n  Outlier information written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def joinTables(self, options):
        """Join tables command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - join_tables] Joining tables containing bin information.')
        self.logger.info('*******************************************************************************')
        
        # read all tables
        headers = {}
        rows = defaultdict(dict)
        binIds = set()
        for f in options.tables:        
            with open(f) as fin:
                headers[f] = [x.strip() for x in fin.readline().split('\t')][1:]
                
                for line in fin:             
                    lineSplit = [x.strip() for x in line.split('\t')]
                    
                    binId = lineSplit[0]
                    binIds.add(binId)
                    
                    for i, header in enumerate(headers[f]):
                        rows[binId][header] = lineSplit[i+1]
                           
        # write merge table  
        oldStdOut = reassignStdOut(options.file)
        
        row = 'Bin Id'        
        for f in options.tables:
            row += '\t' + '\t'.join(headers[f])
        print(row)
              
        for binId in binIds:
            row = binId
            for f in options.tables:
                for header in headers[f]:     
                    row += '\t' + rows[binId].get(header, '')
            print(row)
   
        restoreStdOut(options.file, oldStdOut) 
            
        if options.file:    
            self.logger.info('\n  Joined table written to: ' + options.file)
               
        self.timeKeeper.printTimeStamp()
        
    def modify(self, options):
        """Modify command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - modify] Modifying sequences in bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        if not (options.add or options.remove or options.outlier_file):
            self.logger.warning('  [Warning] No modification to bin requested.\n')
            sys.exit()
            
        if (options.add or options.remove) and options.outlier_file:
            self.logger.warning("  [Warning] The 'outlier_file' option cannot be specified with 'add' or 'remove'.\n")
            sys.exit()
                        
        binTools = BinTools()
        
        if options.add or options.remove:
            binTools.modify(options.bin_file, options.seq_file, options.add, options.remove, options.output_file)
        elif options.outlier_file:
            binTools.removeOutliers(options.bin_file, options.outlier_file, options.output_file)
        
        self.logger.info('  Modified bin written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def unique(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - unique] Ensuring no sequences are assigned to multiple bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
  
        binFiles = self.binFiles(options.bin_folder, options.extension)
        
        binTools = BinTools()
        binTools.unique(binFiles)
            
        self.timeKeeper.printTimeStamp()
        
    def ssuFinder(self, options):
        """SSU finder command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - ssu_finder] Identifying SSU (16S/18S) rRNAs in sequences.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
  
        binFiles = self.binFiles(options.bin_folder, options.extension) 
        
        checkFileExists(options.seq_file)
        makeSurePathExists(options.out_folder)
        
        ssuFinder = SSU_Finder(options.threads)
        ssuFinder.run(options.seq_file, binFiles, options.out_folder, options.evalue, options.concatenate)
            
        self.timeKeeper.printTimeStamp()
        
    def binCompare(self, options):
        """Bin compare command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - bin_compare] Comparing two sets of bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
  
        makeSurePathExists(options.bin_folder1)
        makeSurePathExists(options.bin_folder2)
  
        binFiles1 = self.binFiles(options.bin_folder1, options.extension1)
        binFiles2 = self.binFiles(options.bin_folder2, options.extension2) 
        
        binComparer = BinComparer()
        binComparer.report(binFiles1, binFiles2, options.seq_file, options.output_file)
        
        self.logger.info('')
        self.logger.info('  Detailed bin comparison written to: ' + options.output_file)
       
        self.timeKeeper.printTimeStamp()
        
    def test(self, options):
        """Bin compare command"""
        print ''
        print '*******************************************************************************'
        print '[CheckM - test] Processing E.coli K12-W3310 and verify operation of CheckM.'
        print '*******************************************************************************'
  
        ecoliFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'test_data', '637000110.fna')
        makeSurePathExists(ecoliFile)
        
        testing = Testing()
        
        print '  [Step 1]: Verifying tree command.'
        options.bin_folder = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'test_data')
        options.extension = 'fna'
        options.bQuiet = True
        self.tree(options)
        testing.verifyTree(options.out_folder)
        print '    Passed.'
        
        print '  [Step 2]: Verifying tree_qa command.'
        options.tree_folder = options.out_folder
        options.out_format = 1
        options.file = os.path.join(options.out_folder, 'tree_qa_test.tsv')
        options.bTabTable = True
        self.treeQA(options)
        testing.verifyTreeQA(options.file)
        print '    Passed.'
        
        print '  [Step 3]: Verifying lineage_set command.'
        options.marker_file = os.path.join(options.out_folder, 'lineage_set_test.tsv')
        options.bootstrap = 0
        options.num_genomes_markers = 30
        options.num_genomes_refine = 5
        options.bNoLineageSpecificRefinement = False
        
        options.bRequireTaxonomy = False
        self.lineageSet(options)
        testing.verifyLineageSet(options.marker_file, options.bRequireTaxonomy)
        
        options.bRequireTaxonomy = True
        self.lineageSet(options)
        testing.verifyLineageSet(options.marker_file, options.bRequireTaxonomy)
        print '    Passed.'
        
        print '  [Step 4]: Verifying analyze command.'
        self.analyze(options)
        testing.verifyAnalyze(options.out_folder)
        print '    Passed.'
        
        print '  [Step 5]: Verifying qa command.'
        options.analyze_folder = options.out_folder
        options.out_format = 1
        options.file = os.path.join(options.out_folder, 'qa_test.tsv')
        options.bIndividualMarkers = False
        options.bSkipOrfCorrection = False
        options.bIgnoreThresholds = False
        options.aai_strain = 0.95
        options.e_value = 1e-10
        options.length = 0.7
        options.coverage_file = None
        options.bTabTable = True
        self.qa(options)
        testing.verifyQA(options.file)
        print '    Passed.'

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""    
        self.timeKeeper = TimeKeeper()
         
        try:
            if options.bQuiet:
                logging.basicConfig(format='', level=logging.ERROR)
            else:
                logging.basicConfig(format='', level=logging.INFO)
        except:
            logging.basicConfig(format='', level=logging.INFO)
        
        try:
            if options.file == "stdout":
                options.file = ''
        except:
            pass

        if(options.subparser_name == 'tree'):
            self.tree(options)
        elif(options.subparser_name == 'tree_qa'):
            self.treeQA(options)
        elif(options.subparser_name == 'lineage_set'):
            self.lineageSet(options)
        elif(options.subparser_name == 'taxon_list'):
            self.taxonList(options)
        elif(options.subparser_name == 'taxon_set'):
            self.taxonSet(options)
        elif(options.subparser_name == 'analyze'):
            self.analyze(options)
        elif(options.subparser_name == 'qa'):
            self.qa(options)
        elif(options.subparser_name == 'gc_plot'):
            self.gcPlot(options)
        elif(options.subparser_name == 'coding_plot'):
            self.codingDensityPlot(options)
        elif(options.subparser_name == 'tetra_plot'):
            self.tetraDistPlot(options)
        elif(options.subparser_name == 'dist_plot'):
            self.distributionPlots(options)
        elif(options.subparser_name == 'nx_plot'):
            self.nxPlot(options)
        elif(options.subparser_name == 'len_plot'):
            self.cumulativeLengthPlot(options)
        elif(options.subparser_name == 'len_hist'):
            self.lengthHistogram(options)
        elif(options.subparser_name == 'marker_plot'):
            self.markerPlot(options)
        elif(options.subparser_name == 'par_plot'):
            self.parallelCoordPlot(options)
        elif(options.subparser_name == 'tetra_pca'):
            self.tetraPcaPlot(options)
        elif(options.subparser_name == 'cov_pca'):
            self.coveragePcaPlot(options)
        elif(options.subparser_name == 'gc_bias_plot'):
            self.gcBiasPlot(options)
        elif(options.subparser_name == 'bin_qa_plot'):
            self.binQAPlot(options)
        elif(options.subparser_name == 'unbinned'):
            self.unbinned(options)
        elif(options.subparser_name == 'coverage'):
            self.coverage(options)
        elif(options.subparser_name == 'tetra'):
            self.tetraSignatures(options)
        elif(options.subparser_name == 'profile'):
            self.profile(options)
        elif(options.subparser_name == 'merge'):
            self.merge(options)
        elif(options.subparser_name == 'outliers'):
            self.outliers(options)
        elif(options.subparser_name == 'join_tables'):
            self.joinTables(options)
        elif(options.subparser_name == 'modify'):
            self.modify(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'ssu_finder'):
            self.ssuFinder(options)
        elif(options.subparser_name == 'bin_compare'):
            self.binCompare(options)
        elif(options.subparser_name == 'test'):
            self.test(options)
        else:
            self.logger.error('  [Error] Unknown CheckM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
