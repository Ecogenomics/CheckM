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
import uuid

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
from PCA import PCA
import database as chmdb
from common import makeSurePathExists, checkFileExists, binIdFromFilename

from plot.gcPlots import GcPlots
from plot.codingDensityPlots import CodingDensityPlots
from plot.tetraDistPlots import TetraDistPlots
from plot.distributionPlots import DistributionPlots
from plot.nxPlot import NxPlot
from plot.seqLenPlot import SeqLenPlot
from plot.markerGenePosPlot import MarkerGenePosPlot
from plot.parallelCoordPlot import ParallelCoordPlot
from plot.pcaPlot import PcaPlot

def connectToDatabase(database_name):
    """ Return a database object based on name """
    if database_name is None:
        database_name = os.getenv("CHECKM_DB")
        if database_name is None:
            raise RuntimeError('Cannot connect to DB')

    return chmdb.MarkerDB(db=database_name)

class OptionsParser():
    def __init__(self): pass
    
    def binFiles(self, options):
        targetFiles = []
        if options.bin_folder is not None:
            all_files = os.listdir(options.bin_folder)
            for f in all_files:
                if f.endswith(options.extension):
                    targetFiles.append(os.path.join(options.bin_folder, f))
                    
        if not targetFiles:
            raise "No input files!"
        
        return sorted(targetFiles)

    def analyze(self, options, db=None):
        """Analyze command"""    
        
        targetFiles = self.binFiles(options)    

        # find marker genes in genome bins
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - analyze] Identifying marker genes in bins.'
            print '*******************************************************************************'
        mgf = MarkerGeneFinder(threads=options.threads)
        mgf.find(targetFiles, options.out_folder, options.hmm, bQuiet=options.bQuiet)
        
        self.timeKeeper.printTimeStamp()
        
        # calculate statistics for each genome bin
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - analyze] Calculating genome statistics (e.g., GC, coding density).'
            print '*******************************************************************************'
        binStats = BinStatistics(threads=options.threads)
        binStats.calculate(targetFiles, options.out_folder, bQuiet=options.bQuiet)
        
        self.timeKeeper.printTimeStamp()

    def qa(self, options):
        """QA command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - qa] Tabulating genome statistics.'
            print '*******************************************************************************'

        RP = ResultsParser()
        RP.analyseResults(options.out_folder,
                          options.hmm,
                          eCO=options.e_value,
                          lengthCO=options.length,
                          bQuiet=options.bQuiet,
                          outFile=options.file
                          )
        RP.printSummary(outputFormat=options.out_format, outFile=options.file, bQuiet=options.bQuiet)
                
        if not options.bQuiet and options.file != '':
            print '  QA information written to: ' + options.file
        
        self.timeKeeper.printTimeStamp()

    def align(self, options):
        """Align command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - align] Creating alignments for identified marker genes.'
            print '*******************************************************************************'
                
        if hasattr(options, 'separate'):
            HA = HMMAligner(options.separate, options.consensus, options.out_format)
        else:
            HA = HMMAligner()
            
        bh = False
        if hasattr(options, 'best_alignment'):
            bh=True
            
        HA.makeAlignments(options.out_folder,
                          options.hmm,
                          eCO=options.e_value,
                          lengthCO=options.length,
                          bQuiet=options.bQuiet,
                          bestHit=bh
                          )
        
        self.timeKeeper.printTimeStamp()
        
    def gcPlot(self, options):
        """GC plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - gc_plot] Create GC histogram and delta-GC plot.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plots = GcPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            plots.plot(f, options.distributions)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.gc_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def codingDensityPlot(self, options):
        """Coding density plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - coding_plot] Create coding density (CD) histogram and delta-CD plot.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plots = CodingDensityPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting coding density plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            plots.plot(f, options.distributions)
            
            binId = binIdFromFilename(f) 
            outputFile = os.path.join(options.plot_folder, binId) + '.coding_density_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def tetraDistPlot(self, options):
        """Tetranucleotide distance plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - tetra_plot] Create tetra-distance (TD) histogram and delta-TD plot.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        genomicSignatures = GenomicSignatures()
        tetraSigs = genomicSignatures.read(options.tetra_profile)
            
        plots = TetraDistPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting tetranuclotide distance plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
                
            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.tetra_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def distributionPlots(self, options):
        """Reference distribution plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - dist_plot] Plot GC, CD, and TD distribution plots.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        genomicSignatures = GenomicSignatures()
        tetraSigs = genomicSignatures.read(options.tetra_profile)
            
        plots = DistributionPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting reference distribution plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
                
            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.ref_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def tetraPcaPlot(self, options):
        """PCA plot of tetranucleotide signatures"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - tetra_pca] PCA plot of tetranucleotide signatures.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        print '  Computing PCA of tetranuclotide signatures.\n'
        pca = PCA()
        seqIds, pc, variance = pca.pcaFile(options.tetra_profile, fraction=1.0, bCenter=True, bScale=False)
            
        plots = PcaPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting PCA of tetranuclotide signatures for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
                    
            plots.plot(f, seqIds, pc, variance)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.tetra_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def coveragePcaPlot(self, options):
        """PCA plot of coverage profiles"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - cov_pca] PCA plot of coverage profiles.'
            print '*******************************************************************************'
            
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
            print '  [Warning] Coverage profile is 1 dimensional. PCA requires at least 2 dimensions.'
            sys.exit()
        
        print '  Computing PCA of coverage profiles.\n'
        pca = PCA()
        pc, variance = pca.pcaMatrix(coverageProfiles, fraction=1.0, bCenter=True, bScale=False)
            
        plots = PcaPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting PCA of coverage profiles for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
                    
            plots.plot(f, seqIds, pc, variance)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.cov_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
            
    def nxPlot(self, options):
        """Nx-plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - nx_plot] Creating Nx-plots.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        nx = NxPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting Nx-plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            nx.plot(f)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.nx_plot.' + options.image_type
            nx.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def lenPlot(self, options):
        """Cumulative sequence length plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - len_plot] Creating cumulative sequence length plot.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
            
        plot = SeqLenPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting cumulative sequence length plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            plot.plot(f)
            
            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.seq_len_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def markerPlot(self, options):
        """Marker gene position plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - marker_plot] Creating marker gene position plot.'
            print '*******************************************************************************'
                        
        # generate plot for each bin    
        targetFiles = self.binFiles(options) 
        
        makeSurePathExists(options.plot_folder)
        
        resultsParser = ResultsParser()
        markerGeneStats = resultsParser.parseMarkerGeneStats(options.out_folder)
        binStats = resultsParser.parseBinStatsExt(options.out_folder)
            
        plot = MarkerGenePosPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting marker gene position plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            
            binId = binIdFromFilename(f)
            plot.plot(f, markerGeneStats[binId], binStats[binId])
            
            outputFile = os.path.join(options.plot_folder, binId) + '.marker_pos_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()

    def parallelCoordPlot(self, options):
        """Parallel coordinate plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - par_plot] Parallel coordinate plot of GC and coverage.'
            print '*******************************************************************************'  

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
            if not options.bQuiet:
                print '  Plotting marker gene position plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
                
            binId = binIdFromFilename(f)
            plot.plot(binId, seqStats, coverageStats)
            
            outputFile = os.path.join(options.plot_folder, binId) + '.paralel_coord_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def unbinned(self, options):
        """Unbinned Command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - unbinned] Identify unbinned sequences.'
            print '*******************************************************************************'
  
        binFiles = self.binFiles(options) 
  
        unbinned = Unbinned()
        unbinned.run(binFiles, options.seq_file, options.file, options.min_seq_len, options.bQuiet)
        
        if not options.bQuiet and options.file != '':
            print '  Unbinned sequences written to: ' + options.file
            
        self.timeKeeper.printTimeStamp()
        
    def coverage(self, options):
        """Coverage command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - coverage] Calculate coverage of sequences.'
            print '*******************************************************************************'
    
        binFiles = self.binFiles(options) 
        
        coverage = Coverage(options.threads)
        coverage.calculate(binFiles, options.bam_files, options.output_file, bPairsOnly = options.pairs_only, bQuiet = options.bQuiet)
        
        if not options.bQuiet:
            print '  Coverage information written to: ' + options.output_file    
            
        self.timeKeeper.printTimeStamp()
        
    def tetraSignatures(self, options):
        """Tetranucleotide signature command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - tetra] Calculate tetranucleotide signature of sequences.'
            print '*******************************************************************************'
  
        checkFileExists(options.seq_file)
        
        tetraSig = GenomicSignatures(4, options.threads)
        tetraSig.calculate(options.seq_file, options.output_file, bQuiet = options.bQuiet)
        
        if not options.bQuiet:
            print '  Tetranucletoide signatures written to: ' + options.output_file
            
        self.timeKeeper.printTimeStamp()
        
    def profile(self, options):
        """Profile command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - profile] Calculate percentage of reads mapped to each bin.'
            print '*******************************************************************************'
  
  
        profile = Profile()
        profile.run(options.coverage_file, options.file, bQuiet = options.bQuiet)
        
        if not options.bQuiet and options.file != '':
            print '  Profile information written to: ' + options.file
          
        self.timeKeeper.printTimeStamp()
            
    def outliers(self, options):
        """Outlier command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - outlier] Identify outliers in bins.'
            print '*******************************************************************************'
        
        binFiles = self.binFiles(options) 
        
        binTools = BinTools()
        binTools.identifyOutliers(binFiles, options.distribution, options.report_type, options.output_file, options.bQuiet)

        if not options.bQuiet:
            print '\n  Outlier information written to: ' + options.output_file
            
        self.timeKeeper.printTimeStamp()
        
    def modify(self, options):
        """Modify command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - modify] Modify sequences in a bin.'
            print '*******************************************************************************'
        
        if not (options.add or options.remove or options.outlier_file):
            sys.stderr.write('[Warning] No modification to bin requested.\n')
            sys.exit()
            
        if (options.add or options.remove) and options.outlier_file:
            sys.stderr.write("[Warning] The 'outlier_file' option cannot be specified with 'add' or 'remove'.\n")
            sys.exit()
                        
        binTools = BinTools()
        
        if options.add or options.remove:
            binTools.modify(options.bin_file, options.seq_file, options.add, options.remove, options.output_file)
        elif options.outlier_file:
            binTools.removeOutliers(options.bin_file, options.outlier_file, options.output_file)
        
        if not options.bQuiet:
            print '  Modified bin written to: ' + options.output_file
            
        self.timeKeeper.printTimeStamp()
        
    def unique(self, options):
        """Unique command"""
        print ''
        print '*******************************************************************************'
        print ' [CheckM - unique] Ensure no sequences are assigned to multiple bins.'
        print '*******************************************************************************'
  
        binFiles = self.binFiles(options) 
        
        binTools = BinTools()
        binTools.unique(binFiles)
            
        self.timeKeeper.printTimeStamp()
        
    def makeDB(self, options):
        DB = chmdb.MarkerDB()
        DB.makeDB(options)

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""
        try:
            self.timeKeeper = TimeKeeper(options.bQuiet)
        except:
            # no option to make command quiet, so assume we want visible output
            self.timeKeeper = TimeKeeper(False)
        
        try:
            if options.file == "STDOUT":
                options.file = ''
        except:
            pass

        if(options.subparser_name == 'makeDB'):
            self.makeDB(options)
            return 0

        database_name=None
        try:
            if options.subparser_name in ['analyze', 'qa', 'all']: 
                if options.hmm is None and options.taxonomy is not None:
                    if options.database is None:
                        options.database = os.getenv("CHECKM_DB")
                        if options.database is None:
                            raise RuntimeError('Cannot connect to DB')
    
                    database = chmdb.MarkerDB(db=options.database)
                    tmp = os.path.join('/tmp', str(uuid.uuid4()))
                    database.generateModelFiles(options.taxonomy, tmp)
                    options.hmm = tmp
        except AttributeError, e:
            raise e

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
            self.lenPlot(options)
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
        else:
            sys.stderr.write('[Error] Unknown CheckM command: ' + options.subparser_name + '\n')

        if database_name is not None:
            os.remove(tmp)

        return 0
