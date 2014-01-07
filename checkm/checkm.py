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

import numpy as np

from timeKeeper import TimeKeeper
from resultsParser import ResultsParser, HMMAligner
from markerGeneFinder import MarkerGeneFinder
from binStatistics import BinStatistics
from coverage import Coverage
from genomicSignatures import GenomicSignatures
from unbinned import Unbinned
from profile import Profile
from binTools import BinTools
from ssuFinder import SSU_Finder
from PCA import PCA
from common import makeSurePathExists, checkFileExists, binIdFromFilename

from plot.gcPlots import GcPlots
from plot.codingDensityPlots import CodingDensityPlots
from plot.tetraDistPlots import TetraDistPlots
from plot.distributionPlots import DistributionPlots
from plot.nxPlot import NxPlot
from plot.cumulativeLengthPlot import CumulativeLengthPlot
from plot.lengthHistogram import LengthHistogram
from plot.markerGenePosPlot import MarkerGenePosPlot
from plot.parallelCoordPlot import ParallelCoordPlot
from plot.pcaPlot import PcaPlot

class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()   
    
    def binFiles(self, options):
        binFiles = []
        if options.bin_folder is not None:
            all_files = os.listdir(options.bin_folder)
            for f in all_files:
                if f.endswith(options.extension):
                    binFiles.append(os.path.join(options.bin_folder, f))
                    
        if not binFiles:
            self.logger.error("  [Error] No bins found. Check the extension used to identify bins.")
            sys.exit()
        
        return sorted(binFiles)

    def analyze(self, options, db=None):
        """Analyze command"""   
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - analyze] Identifying marker genes in bins.')
        self.logger.info('*******************************************************************************')
                
        binFiles = self.binFiles(options)   
        makeSurePathExists(options.out_folder) 

        # find marker genes in genome bins    
        mgf = MarkerGeneFinder(options.threads)
        mgf.find(binFiles, options.out_folder, options.marker_file)
        
        self.timeKeeper.printTimeStamp()
        
        # calculate statistics for each genome bin
        self.logger.info('')
        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.out_folder)
        
        self.timeKeeper.printTimeStamp()

    def qa(self, options):
        """QA command"""
        if options.file != '':
            self.logger.info('')
            self.logger.info('*******************************************************************************')
            self.logger.info(' [CheckM - qa] Tabulating genome statistics.')
            self.logger.info('*******************************************************************************')

        RP = ResultsParser()
        RP.analyseResults(options.out_folder,
                          options.marker_file,
                          bIgnoreThresholds = options.bIgnoreThresholds,
                          evalueThreshold = options.e_value,
                          lengthThreshold = options.length,
                          bSkipOrfCorrection = options.bSkipOrfCorrection
                          )
        RP.printSummary(options.out_format, options.bTabTable, options.file)
                
        if options.file != '':
            print '  QA information written to: ' + options.file

    def align(self, options):
        """Align command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - align] Creating alignments for identified marker genes.')
        self.logger.info('*******************************************************************************')
                
        if hasattr(options, 'separate'):
            HA = HMMAligner(options.separate, options.consensus, options.out_format)
        else:
            HA = HMMAligner()
            
        bh = False
        if hasattr(options, 'best_alignment'):
            bh=True
            
        HA.makeAlignments(options.out_folder,
                          options.hmm,
                          evalueThreshold = options.e_value,
                          lengthThreshold = options.length,
                          bestHit=bh
                          )
        
        self.timeKeeper.printTimeStamp()
        
    def gcPlot(self, options):
        """GC plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - gc_plot] Create GC histogram and delta-GC plot.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plots = GcPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
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
        self.logger.info(' [CheckM - coding_plot] Create coding density (CD) histogram and delta-CD plot.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plots = CodingDensityPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting coding density plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
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
        self.logger.info(' [CheckM - tetra_plot] Create tetra-distance (TD) histogram and delta-TD plot.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)
            
        plots = TetraDistPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting tetranuclotide distance plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
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
        self.logger.info(' [CheckM - dist_plot] Plot GC, CD, and TD distribution plots.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        genomicSignatures = GenomicSignatures()
        tetraSigs = genomicSignatures.read(options.tetra_profile)
            
        plots = DistributionPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting reference distribution plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
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
        self.logger.info(' [CheckM - tetra_pca] PCA plot of tetranucleotide signatures.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        self.logger.info('  Computing PCA of tetranuclotide signatures.\n')
        pca = PCA()
        seqIds, pc, variance = pca.pcaFile(options.tetra_profile, fraction=1.0, bCenter=True, bScale=False)
            
        plots = PcaPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting PCA of tetranuclotide signatures for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
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
        self.logger.info(' [CheckM - cov_pca] PCA plot of coverage profiles.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        coverage = Coverage()
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
        
        print '  Computing PCA of coverage profiles.\n'
        pca = PCA()
        pc, variance = pca.pcaMatrix(coverageProfiles, fraction=1.0, bCenter=True, bScale=False)
            
        plots = PcaPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting PCA of coverage profiles for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
            filesProcessed += 1
                    
            plots.plot(f, seqIds, pc, variance)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.cov_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
            
    def nxPlot(self, options):
        """Nx-plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - nx_plot] Creating Nx-plots.')
        self.logger.info('*******************************************************************************')
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        nx = NxPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting Nx-plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
            filesProcessed += 1
            nx.plot(f)
            
            binId = binIdFromFilename(f)
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
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plot = CumulativeLengthPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting cumulative sequence length plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
            filesProcessed += 1
            plot.plot(f)
            
            binId = binIdFromFilename(f)
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
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plot = LengthHistogram(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting sequence length histogram for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
            filesProcessed += 1
            plot.plot(f)
            
            binId = binIdFromFilename(f)
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
                        
        # generate plot for each bin    
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        resultsParser = ResultsParser()
        markerGeneStats = resultsParser.parseMarkerGeneStats(options.out_folder)
        binStats = resultsParser.parseBinStatsExt(options.out_folder)
            
        plot = MarkerGenePosPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting marker gene position plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
            filesProcessed += 1
            
            binId = binIdFromFilename(f)
            plot.plot(f, markerGeneStats[binId], binStats[binId])
            
            outputFile = os.path.join(options.plot_folder, binId) + '.marker_pos_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()

    def parallelCoordPlot(self, options):
        """Parallel coordinate plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - par_plot] Parallel coordinate plot of GC and coverage.')
        self.logger.info('*******************************************************************************')  

        targetFiles = self.binFiles(options) 

        makeSurePathExists(options.plot_folder)
        checkFileExists(options.coverage_file)

        # read sequence stats file
        resultsParser = ResultsParser()
        seqStats = resultsParser.parseSeqStats(options.out_folder)

        # read coverage stats file
        coverage = Coverage()
        coverageStats = coverage.parseCoverage(options.coverage_file)
        
        # create plot for each bin
        plot = ParallelCoordPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            self.logger.info('  Plotting marker gene position plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles)))
            filesProcessed += 1
                
            binId = binIdFromFilename(f)
            plot.plot(binId, seqStats, coverageStats)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.paralel_coord_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def unbinned(self, options):
        """Unbinned Command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - unbinned] Identify unbinned sequences.')
        self.logger.info('*******************************************************************************')
  
        binFiles = self.binFiles(options) 
  
        unbinned = Unbinned()
        unbinned.run(binFiles, options.seq_file, options.file, options.min_seq_len)
        
        if options.file != '':
            self.logger.info('  Unbinned sequences written to: ' + options.file)
            
        self.timeKeeper.printTimeStamp()
        
    def coverage(self, options):
        """Coverage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - coverage] Calculate coverage of sequences.')
        self.logger.info('*******************************************************************************')
    
        binFiles = self.binFiles(options) 
        
        coverage = Coverage(options.threads)
        coverage.run(binFiles, options.bam_files, options.output_file, options.all_reads, options.min_align, options.max_edit_dist)
        
        self.logger.info('  Coverage information written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def tetraSignatures(self, options):
        """Tetranucleotide signature command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra] Calculate tetranucleotide signature of sequences.')
        self.logger.info('*******************************************************************************')
  
        checkFileExists(options.seq_file)
        
        tetraSig = GenomicSignatures(4, options.threads)
        tetraSig.calculate(options.seq_file, options.output_file)
        
        self.logger.info('  Tetranucletoide signatures written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def profile(self, options):
        """Profile command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - profile] Calculate percentage of reads mapped to each bin.')
        self.logger.info('*******************************************************************************')
  
        profile = Profile()
        profile.run(options.coverage_file, options.file)
        
        if options.file != '':
            self.logger.info('  Profile information written to: ' + options.file)
          
        self.timeKeeper.printTimeStamp()
            
    def outliers(self, options):
        """Outlier command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - outlier] Identify outliers in bins.')
        self.logger.info('*******************************************************************************')
        
        binFiles = self.binFiles(options) 
        
        binTools = BinTools()
        binTools.identifyOutliers(binFiles, options.distribution, options.report_type, options.output_file)

        self.logger.info('\n  Outlier information written to: ' + options.output_file)
            
        self.timeKeeper.printTimeStamp()
        
    def modify(self, options):
        """Modify command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - modify] Modify sequences in a bin.')
        self.logger.info('*******************************************************************************')
        
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
        self.logger.info('[CheckM - unique] Ensure no sequences are assigned to multiple bins.')
        self.logger.info('*******************************************************************************')
  
        binFiles = self.binFiles(options) 
        
        binTools = BinTools()
        binTools.unique(binFiles)
            
        self.timeKeeper.printTimeStamp()
        
    def ssuFinder(self, options):
        """SSU finder command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - ssu_finder] Identify SSU (16S/18S) rRNAs in sequences.')
        self.logger.info('*******************************************************************************')
  
        binFiles = self.binFiles(options) 
        
        checkFileExists(options.seq_file)
        makeSurePathExists(options.out_folder)
        
        ssuFinder = SSU_Finder(options.threads)
        ssuFinder.run(options.seq_file, binFiles, options.out_folder, options.evalue, options.concatenate)
            
        self.timeKeeper.printTimeStamp()

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""    
        self.timeKeeper = TimeKeeper()
         
        if options.bQuiet:
            logging.basicConfig(format='', level=logging.ERROR)
        else:
            logging.basicConfig(format='', level=logging.INFO)
        
        try:
            if options.file == "STDOUT":
                options.file = ''
        except:
            pass

        if(options.subparser_name == 'analyze'):
            self.analyze(options)
        elif(options.subparser_name == 'qa'):
            self.qa(options)
        elif(options.subparser_name == 'all'):
            self.analyze(options)
            self.qa(options)
            self.align(options)
        elif(options.subparser_name == 'align'):
            self.align(options)
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
        elif(options.subparser_name == 'unbinned'):
            self.unbinned(options)
        elif(options.subparser_name == 'coverage'):
            self.coverage(options)
        elif(options.subparser_name == 'tetra'):
            self.tetraSignatures(options)
        elif(options.subparser_name == 'profile'):
            self.profile(options)
        elif(options.subparser_name == 'outliers'):
            self.outliers(options)
        elif(options.subparser_name == 'modify'):
            self.modify(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'ssu_finder'):
            self.ssuFinder(options)
        else:
            self.logger.error('  [Error] Unknown CheckM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
