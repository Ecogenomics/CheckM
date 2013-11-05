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
import uuid

from timeKeeper import TimeKeeper
from resultsParser import ResultsParser, HMMAligner
from markerGeneFinder import MarkerGeneFinder
from binStatistics import BinStatistics
from coverage import Coverage
from unbinned import Unbinned
from profile import Profile
import database as chmdb
import defaultValues

from plot.gcPlots import gcPlots
from plot.nxPlot import nxPlot
from plot.seqLenPlot import seqLenPlot
from plot.markerGenePosPlot import markerGenePosPlot

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
        
        return targetFiles

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
        
        # save marker gene stats needs for down-stream plots
        RP.printSummary(outputFormat=8, outFile=options.out_folder + '/' + defaultValues.__CHECKM_DEFAULT_MARKER_GENE_STATS__, bQuiet=True)
        
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
        
    def gc_plot(self, options):
        """GC plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - gc_plot] Creating GC histograms.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        if not os.path.exists(options.plot_folder):
            os.mkdir(options.plot_folder)
            
        plots = gcPlots(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            plots.plot(f)
            
            outputFile = os.path.join(options.plot_folder, os.path.basename(f[0:f.rfind('.')])) + '.gc_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
            
    def nx_plot(self, options):
        """Nx-plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - nx_plot] Creating Nx-plots.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        if not os.path.exists(options.plot_folder):
            os.mkdir(options.plot_folder)
            
        nx = nxPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting Nx-plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            nx.plot(f)
            
            outputFile = os.path.join(options.plot_folder, os.path.basename(f[0:f.rfind('.')])) + '.nx_plot.' + options.image_type
            nx.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def len_plot(self, options):
        """Cumulative sequence length plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - len_plot] Creating cumulative sequence length plot.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        if not os.path.exists(options.plot_folder):
            os.mkdir(options.plot_folder)
            
        plot = seqLenPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting cumulative sequence length plot for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            plot.plot(f)
            
            outputFile = os.path.join(options.plot_folder, os.path.basename(f[0:f.rfind('.')])) + '.seq_len_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            if not options.bQuiet:
                print '    Plot written to: ' + outputFile
            
        self.timeKeeper.printTimeStamp()
        
    def marker_plot(self, options):
        """Marker gene position plot command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - marker_plot] Creating marker gene position plot.'
            print '*******************************************************************************'
            
        targetFiles = self.binFiles(options) 
        
        if not os.path.exists(options.plot_folder):
            os.mkdir(options.plot_folder)
            
        plot = markerGenePosPlot(options)
        filesProcessed = 1
        for f in targetFiles:  
            if not options.bQuiet:
                print '  Plotting marker gene positions for %s (%d of %d)' % (f, filesProcessed, len(targetFiles))
                filesProcessed += 1
            plot.plot(f, options.results_folder)
            
            outputFile = os.path.join(options.plot_folder, os.path.basename(f[0:f.rfind('.')])) + '.marker_position_plot.' + options.image_type
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
        """Coverage Command"""
        if not options.bQuiet:
            print ''
            print '*******************************************************************************'
            print ' [CheckM - coverage] Calculate coverage of sequences.'
            print '*******************************************************************************'
  
        
        binFiles = self.binFiles(options) 
        
        coverage = Coverage(threads = options.threads)
        coverage.calculate(binFiles, options.bam_files, options.file, bPairsOnly = options.pairs_only, bQuiet = options.bQuiet)
        
        if not options.bQuiet and options.file != '':
            print '  Coverage information written to: ' + options.file
            
        self.timeKeeper.printTimeStamp()
        
    def profile(self, options):
        """Profile Command"""
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
        
    def makeDB(self, options):
        DB = chmdb.MarkerDB()
        DB.makeDB(options)

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""
        self.timeKeeper = TimeKeeper(options.bQuiet)
        
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
            self.gc_plot(options)
        elif(options.subparser_name == 'nx_plot'):
            self.nx_plot(options)
        elif(options.subparser_name == 'len_plot'):
            self.len_plot(options)
        elif(options.subparser_name == 'marker_plot'):
            self.marker_plot(options)
        elif(options.subparser_name == 'unbinned'):
            self.unbinned(options)
        elif(options.subparser_name == 'coverage'):
            self.coverage(options)
        elif(options.subparser_name == 'profile'):
            self.profile(options)

        if database_name is not None:
            os.remove(tmp)

        return 0
