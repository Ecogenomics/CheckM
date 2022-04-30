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

from checkm.defaultValues import DefaultValues
from checkm.timeKeeper import TimeKeeper
from checkm.markerSets import MarkerSetParser, BinMarkerSets
from checkm.resultsParser import ResultsParser
from checkm.hmmerAligner import HmmerAligner
from checkm.markerGeneFinder import MarkerGeneFinder
from checkm.pplacer import PplacerRunner
from checkm.treeParser import TreeParser
from checkm.taxonParser import TaxonParser
from checkm.aminoAcidIdentity import AminoAcidIdentity
from checkm.binStatistics import BinStatistics
from checkm.coverage import Coverage
from checkm.coverageWindows import CoverageWindows
from checkm.genomicSignatures import GenomicSignatures
from checkm.unbinned import Unbinned
from checkm.merger import Merger
from checkm.profile import Profile
from checkm.binTools import BinTools
from checkm.ssuFinder import SSU_Finder
from checkm.common import (checkBinInputExists,
                           makeSurePathExists,
                           checkFileExists,
                           binIdFromFilename,
                           getBinIdsFromOutDir,
                           checkDirExists)

from checkm.plot.gcPlots import GcPlots
from checkm.plot.codingDensityPlots import CodingDensityPlots
from checkm.plot.tetraDistPlots import TetraDistPlots
from checkm.plot.distributionPlots import DistributionPlots
from checkm.plot.nxPlot import NxPlot
from checkm.plot.lengthHistogram import LengthHistogram
from checkm.plot.markerGenePosPlot import MarkerGenePosPlot
from checkm.plot.gcBiasPlots import GcBiasPlot

from checkm.util.seqUtils import checkNuclotideSeqs, checkProteinSeqs

from checkm.checkmData import DBManager

from checkm.test.test_ecoli import VerifyEcoli


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger('timestamp')
        self.timeKeeper = TimeKeeper()

    def updateCheckM_DB(self, options):
        self.logger.info(
            '[CheckM - data] Check for database updates. [%s]' % options.action[0])

        action = options.action
        if action and action[0] == 'setRoot':
            if len(action) > 2:
                DBM = DBManager(set_path=action[1], configFile=action[2])
            elif len(action) > 1:
                DBM = DBManager(set_path=action[1])
        else:
            self.logger.error(
                'Path to the CheckM reference data must be specified.')

    def binFiles(self, binInput, binExtension, bCalledGenes):
        binFiles = []
        binIDs = set()
        isInputDir = True
        if binInput is not None:
            if os.path.isdir(binInput):
                if binExtension[0] != '.':
                    binExtension = '.' + binExtension

                all_files = os.listdir(binInput)
                for f in all_files:
                    if f.endswith(binExtension):
                        binFile = os.path.join(binInput, f)
                        if os.stat(binFile).st_size == 0:
                            self.logger.warning(
                                "Skipping bin %s as it has a size of 0 bytes." % f)
                        else:
                            binFiles.append(binFile)
                            binIDs.add(os.path.basename(binFile))
            else:
                with open(binInput, "r") as oh:
                    for line in oh:
                        files = line.strip().split("\t")
                        binFile = files[1]
                        if bCalledGenes:
                            binFile = files[2]
                        if not os.path.exists(binFile):
                            self.logger.warning(
                                "Skipping bin %s as it doesn't exists." % binFile)
                        elif os.stat(binFile).st_size == 0:
                            self.logger.warning(
                                "Skipping bin %s as it has a size of 0 bytes." % binFile)
                        else:
                            binFiles.append(binFile)
                            binIDs.add(os.path.basename(binFile))

        if not binFiles:
            if isInputDir:
                self.logger.error(
                    "No bins found. Check the extension (-x) used to identify bins.")
            else:
                self.logger.error(
                    "No bins found. Check the bins input table to verify bins exists.")
            sys.exit(1)

        if len(binIDs) != len(binFiles):
            self.logger.error(
                "There are redundant bin IDs, please check and update.")
            sys.exit(1)

        return sorted(binFiles)

    def tree(self, options):
        """Tree command"""

        self.logger.info(
            '[CheckM - tree] Placing bins in reference genome tree.')

        binFiles = self.binFiles(
            options.bin_input, options.extension, options.bCalledGenes)

        if not options.bCalledGenes:
            if not checkNuclotideSeqs(binFiles):
                return
        else:
            if not checkProteinSeqs(binFiles):
                return

        # setup directory structure
        makeSurePathExists(os.path.join(options.output_dir, 'bins'))
        makeSurePathExists(os.path.join(options.output_dir, 'storage'))

        # find phylogenetically informative genes in genome bins
        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(binFiles,
                                 options.output_dir,
                                 DefaultValues.HMMER_TABLE_PHYLO_OUT,
                                 DefaultValues.HMMER_PHYLO_OUT,
                                 DefaultValues.PHYLO_HMM_MODELS,
                                 options.bKeepAlignment,
                                 options.bNucORFs,
                                 options.bCalledGenes)

        # write model information to file
        markerSetParser = MarkerSetParser(options.threads)
        hmmModelInfoFile = os.path.join(
            options.output_dir, 'storage', DefaultValues.PHYLO_HMM_MODEL_INFO)
        markerSetParser.writeBinModels(binIdToModels, hmmModelInfoFile)

        # calculate statistics for each genome bin

        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.output_dir,
                           DefaultValues.BIN_STATS_PHYLO_OUT)

        # align identified marker genes

        HA = HmmerAligner(options.threads)
        resultsParser = HA.makeAlignmentToPhyloMarkers(options.output_dir,
                                                       DefaultValues.PHYLO_HMM_MODELS,
                                                       DefaultValues.HMMER_TABLE_PHYLO_OUT,
                                                       binIdToModels,
                                                       False,
                                                       DefaultValues.E_VAL,
                                                       DefaultValues.LENGTH,
                                                       False,
                                                       os.path.join(
                                                           options.output_dir, 'storage', 'tree')
                                                       )

        # place bins into genome tree

        # fix at one thread to keep memory requirements reasonable
        pplacer = PplacerRunner(threads=options.pplacer_threads)
        pplacer.run(binFiles, resultsParser,
                    options.output_dir, options.bReducedTree)

        self.timeKeeper.printTimeStamp()

    def treeQA(self, options):
        """QA command"""

        self.logger.info(
            '[CheckM - tree_qa] Assessing phylogenetic markers found in each bin.')

        checkDirExists(options.tree_dir)

        # set HMM file for each bin
        markerSetParser = MarkerSetParser()
        hmmModelInfoFile = os.path.join(
            options.tree_dir, 'storage', DefaultValues.PHYLO_HMM_MODEL_INFO)
        binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        # calculate marker gene statistics
        RP = ResultsParser(binIdToModels)
        binStats = RP.analyseResults(options.tree_dir,
                                     DefaultValues.BIN_STATS_PHYLO_OUT,
                                     DefaultValues.HMMER_TABLE_PHYLO_OUT)

        # determine taxonomy of each bin

        treeParser = TreeParser()
        treeParser.printSummary(
            options.out_format, options.tree_dir, RP, options.bTabTable, options.file, binStats)

        if options.file != '':
            self.logger.info('QA information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def lineageSet(self, options, db=None):
        """Lineage set command"""

        self.logger.info(
            '[CheckM - lineage_set] Inferring lineage-specific marker sets.')

        checkDirExists(options.tree_dir)

        # set HMM file for each bin
        markerSetParser = MarkerSetParser()
        hmmModelInfoFile = os.path.join(
            options.tree_dir, 'storage', DefaultValues.PHYLO_HMM_MODEL_INFO)
        binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        # calculate marker gene statistics
        resultsParser = ResultsParser(binIdToModels)
        resultsParser.analyseResults(options.tree_dir,
                                     DefaultValues.BIN_STATS_PHYLO_OUT,
                                     DefaultValues.HMMER_TABLE_PHYLO_OUT)

        # These options are incompatible with how the lineage-specific marker set is selected, so
        # the default values are currently hard-coded

        options.num_genomes_markers = 2
        options.bootstrap = 0
        options.bRequireTaxonomy = False

        treeParser = TreeParser()
        treeParser.getBinMarkerSets(options.tree_dir, options.marker_file,
                                    options.num_genomes_markers,
                                    options.bootstrap, options.bNoLineageSpecificRefinement,
                                    options.bForceDomain, options.bRequireTaxonomy,
                                    resultsParser, options.unique, options.multi)

        self.logger.info('Marker set written to: ' + options.marker_file)

        self.timeKeeper.printTimeStamp()

    def taxonList(self, options, db=None):
        """Lineage set command"""

        self.logger.info(
            '[CheckM - taxon_list] Listing available taxonomic-specific marker sets.')

        taxonParser = TaxonParser()
        taxonParser.list(options.rank)

        self.timeKeeper.printTimeStamp()

    def taxonSet(self, options, db=None):
        """Taxon set command"""

        self.logger.info(
            '[CheckM - taxon_set] Generate taxonomic-specific marker set.')

        path = os.path.split(options.marker_file)[0]
        if path:
            makeSurePathExists(path)

        taxonParser = TaxonParser()
        bValidSet = taxonParser.markerSet(
            options.rank, options.taxon, options.marker_file)

        if bValidSet:
            self.logger.info('Marker set written to: ' + options.marker_file)

        self.timeKeeper.printTimeStamp()

    def analyze(self, options, db=None):
        """Analyze command"""

        self.logger.info(
            '[CheckM - analyze] Identifying marker genes in bins.')

        binFiles = self.binFiles(
            options.bin_input, options.extension, options.bCalledGenes)

        if not options.bCalledGenes:
            if not checkNuclotideSeqs(binFiles):
                return
        else:
            if not checkProteinSeqs(binFiles):
                return

        # setup directory structure
        makeSurePathExists(options.output_dir)
        makeSurePathExists(os.path.join(options.output_dir, 'bins'))
        makeSurePathExists(os.path.join(options.output_dir, 'storage'))
        makeSurePathExists(os.path.join(
            options.output_dir, 'storage', 'aai_qa'))

        # find marker genes in genome bins
        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(binFiles,
                                 options.output_dir,
                                 DefaultValues.HMMER_TABLE_OUT,
                                 DefaultValues.HMMER_OUT,
                                 options.marker_file,
                                 options.bKeepAlignment,
                                 options.bNucORFs,
                                 options.bCalledGenes)

        markerSetParser = MarkerSetParser(options.threads)
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.output_dir,
                                                             getBinIdsFromOutDir(
                                                                 options.output_dir),
                                                             options.marker_file)

        hmmModelInfoFile = os.path.join(
            options.output_dir, 'storage', DefaultValues.CHECKM_HMM_MODEL_INFO)
        markerSetParser.writeBinModels(binIdToModels, hmmModelInfoFile)

        self.timeKeeper.printTimeStamp()

        # HMM model file
        if markerSetParser.markerFileType(options.marker_file) == BinMarkerSets.HMM_MODELS_SET:
            markerFile = options.marker_file
        else:
            markerFile = DefaultValues.HMM_MODELS

        # align marker genes with multiple hits within a bin
        HA = HmmerAligner(options.threads)
        HA.makeAlignmentsOfMultipleHits(options.output_dir,
                                        markerFile,
                                        DefaultValues.HMMER_TABLE_OUT,
                                        binIdToModels,
                                        binIdToBinMarkerSets,
                                        False,
                                        DefaultValues.E_VAL,
                                        DefaultValues.LENGTH,
                                        os.path.join(
                                            options.output_dir, 'storage', 'aai_qa')
                                        )

        self.timeKeeper.printTimeStamp()

        # calculate statistics for each genome bin
        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.output_dir,
                           DefaultValues.BIN_STATS_OUT)

        self.timeKeeper.printTimeStamp()

        # align top hit to each marker if requested
        if options.bAlignTopHit:
            alignmentOutputFolder = os.path.join(
                options.output_dir, 'storage', 'alignments')
            makeSurePathExists(alignmentOutputFolder)

            HA = HmmerAligner(options.threads)
            resultsParser = HA.makeAlignmentTopHit(options.output_dir,
                                                   options.marker_file,
                                                   DefaultValues.HMMER_TABLE_OUT,
                                                   binIdToModels,
                                                   False,
                                                   DefaultValues.E_VAL,
                                                   DefaultValues.LENGTH,
                                                   True,
                                                   alignmentOutputFolder
                                                   )

            # report marker gene data
            fout = open(os.path.join(alignmentOutputFolder,
                        'alignment_info.tsv'), 'w')
            fout.write('Marker Id\tLength (bp)\n')
            markerIds = resultsParser.models[list(
                resultsParser.models.keys())[0]].keys()
            for markerId in markerIds:
                fout.write('%s\t%d\n' % (markerId, resultsParser.models[list(
                    resultsParser.models.keys())[0]][markerId].leng))
            fout.close()

            self.logger.info(
                'Alignments to top hits stored in: ' + alignmentOutputFolder)

            self.timeKeeper.printTimeStamp()

    def qa(self, options):
        """QA command"""
        self.logger.info('[CheckM - qa] Tabulating genome statistics.')

        checkDirExists(options.analyze_dir)

        if options.exclude_markers:
            checkFileExists(options.exclude_markers)

        # calculate AAI between marks with multiple hits in a single bin
        aai = AminoAcidIdentity()
        aai.run(options.aai_strain, options.analyze_dir, options.alignment_file)

        # get HMM file for each bin
        markerSetParser = MarkerSetParser(options.threads)

        hmmModelInfoFile = os.path.join(
            options.analyze_dir, 'storage', DefaultValues.CHECKM_HMM_MODEL_INFO)
        binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.analyze_dir,
                                                             getBinIdsFromOutDir(
                                                                 options.analyze_dir),
                                                             options.marker_file,
                                                             options.exclude_markers)

        # get results for each bin
        RP = ResultsParser(binIdToModels)
        RP.analyseResults(options.analyze_dir,
                          DefaultValues.BIN_STATS_OUT,
                          DefaultValues.HMMER_TABLE_OUT,
                          bIgnoreThresholds=options.bIgnoreThresholds,
                          evalueThreshold=options.e_value,
                          lengthThreshold=options.length,
                          bSkipPseudoGeneCorrection=options.bSkipPseudoGeneCorrection,
                          bSkipAdjCorrection=options.bSkipAdjCorrection
                          )

        RP.printSummary(options.out_format,
                        aai, binIdToBinMarkerSets,
                        options.bIndividualMarkers,
                        options.coverage_file,
                        options.bTabTable,
                        options.file,
                        anaFolder=options.analyze_dir)
        RP.cacheResults(options.analyze_dir,
                        binIdToBinMarkerSets,
                        options.bIndividualMarkers)

        if options.file != '':
            self.logger.info('QA information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def gcPlot(self, options):
        """GC plot command"""
        self.logger.info(
            '[CheckM - gc_plot] Creating GC histogram and delta-GC plot.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input,
            options.extension, False)

        plots = GcPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('Plotting GC plots for %s (%d of %d)' %
                             (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, options.distributions)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(
                options.output_dir, binId) + '.gc_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def codingDensityPlot(self, options):
        """Coding density plot command"""
        self.logger.info(
            '[CheckM - coding_plot] Creating coding density histogram and delta-CD plot.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        plots = CodingDensityPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('Plotting coding density plots for %s (%d of %d)' % (
                f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, options.distributions)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(
                options.output_dir, binId) + '.coding_density_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def tetraDistPlot(self, options):
        """Tetranucleotide distance plot command"""
        self.logger.info(
            '[CheckM - tetra_plot] Creating tetra-distance histogram and delta-TD plot.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)

        plots = TetraDistPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('Plotting tetranuclotide distance plots for %s (%d of %d)' % (
                f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)

            outputFile = os.path.join(
                options.output_dir, binId) + '.tetra_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def distributionPlots(self, options):
        """Reference distribution plot command"""
        self.logger.info(
            '[CheckM - dist_plot] Creating GC, CD, and TD distribution plots.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)

        plots = DistributionPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('Plotting reference distribution plots for %s (%d of %d)' % (
                f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)

            outputFile = os.path.join(
                options.output_dir, binId) + '.ref_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def gcBiasPlot(self, options):
        """GC bias plot command"""

        self.logger.info(
            '[CheckM - gc_bias_plot] Plotting bin coverage as a function of GC.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        coverageWindows = CoverageWindows(options.threads)
        coverageProfile = coverageWindows.run(
            binFiles, options.bam_file, options.all_reads, options.min_align, options.max_edit_dist, options.window_size)

        plots = GcBiasPlot(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('Plotting GC plots for %s (%d of %d)' %
                             (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, coverageProfile)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(
                options.output_dir, binId) + '.gc_bias_plot.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def nxPlot(self, options):
        """Nx-plot command"""

        self.logger.info('[CheckM - nx_plot] Creating Nx-plots.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        nx = NxPlot(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('Plotting Nx-plot for %s (%d of %d)' %
                             (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            nx.plot(f)

            outputFile = os.path.join(
                options.output_dir, binId) + '.nx_plot.' + options.image_type
            nx.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def lengthHistogram(self, options):
        """Sequence length histogram command"""

        self.logger.info(
            '[CheckM - len_hist] Creating sequence length histogram.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        plot = LengthHistogram(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('Plotting sequence length histogram for %s (%d of %d)' % (
                binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            plot.plot(f)

            outputFile = os.path.join(
                options.output_dir, binId) + '.len_hist.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def markerPlot(self, options):
        """Marker gene position plot command"""

        self.logger.info(
            '[CheckM - marker_plot] Creating marker gene position plot.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(options.output_dir)

        # generate plot for each bin
        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        resultsParser = ResultsParser(None)
        markerGeneStats = resultsParser.parseMarkerGeneStats(
            options.results_dir)
        binStats = resultsParser.parseBinStatsExt(options.results_dir)

        plot = MarkerGenePosPlot(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('Plotting marker gene position plot for %s (%d of %d)' % (
                binId, filesProcessed, len(binFiles)))
            filesProcessed += 1

            if binId not in markerGeneStats or binId not in binStats:
                continue  # bin has no marker genes

            bPlotted = plot.plot(f, markerGeneStats[binId], binStats[binId])

            if bPlotted:
                outputFile = os.path.join(
                    options.output_dir, binId) + '.marker_pos_plot.' + options.image_type
                plot.savePlot(outputFile, dpi=options.dpi)
                self.logger.info('Plot written to: ' + outputFile)
            else:
                self.logger.info('No marker genes found in bin.')

        self.timeKeeper.printTimeStamp()

    def unbinned(self, options):
        """Unbinned Command"""

        self.logger.info('[CheckM - unbinned] Identify unbinned sequences.')

        checkBinInputExists(options.bin_input, False)

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        unbinned = Unbinned()
        unbinned.run(binFiles, options.seq_file, options.output_seq_file,
                     options.output_stats_file, options.min_seq_len)

        self.logger.info('Unbinned sequences written to: ' +
                         options.output_seq_file)
        self.logger.info(
            'Unbinned sequences statistics written to: ' + options.output_stats_file)

        self.timeKeeper.printTimeStamp()

    def coverage(self, options):
        """Coverage command"""

        self.logger.info(
            '[CheckM - coverage] Calculating coverage of sequences.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(os.path.dirname(options.output_file))

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        coverage = Coverage(options.threads)
        coverage.run(binFiles, options.bam_files, options.output_file, options.all_reads,
                     options.min_align, options.max_edit_dist, options.min_qc)

        self.logger.info(
            'Coverage information written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def tetraSignatures(self, options):
        """Tetranucleotide signature command"""

        self.logger.info(
            '[CheckM - tetra] Calculating tetranucleotide signature of sequences.')

        checkFileExists(options.seq_file)
        makeSurePathExists(os.path.dirname(options.output_file))

        tetraSig = GenomicSignatures(4, options.threads)
        tetraSig.calculate(options.seq_file, options.output_file)

        self.logger.info(
            'Tetranucletoide signatures written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def profile(self, options):
        """Profile command"""

        self.logger.info(
            '[CheckM - profile] Calculating percentage of reads mapped to each bin.')

        checkFileExists(options.coverage_file)

        profile = Profile()
        profile.run(options.coverage_file, options.file, options.bTabTable)

        if options.file != '':
            self.logger.info('Profile information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def merge(self, options):
        """Merge command"""

        self.logger.info(
            '[CheckM - merge] Identifying bins with complementary sets of marker genes.')

        checkBinInputExists(options.bin_input, options.bCalledGenes)

        binFiles = self.binFiles(
            options.bin_input, options.extension, options.bCalledGenes)

        if not options.bCalledGenes:
            if not checkNuclotideSeqs(binFiles):
                return
        else:
            if not checkProteinSeqs(binFiles):
                return

        markerSetParser = MarkerSetParser()
        if markerSetParser.markerFileType(options.marker_file) == BinMarkerSets.TREE_MARKER_SET:
            self.logger.error(
                'Merge command requires a taxonomic-specific marker set or a user-defined HMM file.\n')
            return

        # setup directory structure
        makeSurePathExists(options.output_dir)
        makeSurePathExists(os.path.join(options.output_dir, 'bins'))
        makeSurePathExists(os.path.join(options.output_dir, 'storage'))
        makeSurePathExists(os.path.join(options.output_dir, 'storage', 'hmms'))

        binIds = []
        for binFile in binFiles:
            binIds.append(binIdFromFilename(binFile))

        # find marker genes in genome bins
        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(binFiles,
                                 options.output_dir,
                                 "merger.table.txt",
                                 "merger.hmmer3",
                                 options.marker_file,
                                 False,
                                 False,
                                 options.bCalledGenes)

        # get HMM file for each bin
        markerSetParser = MarkerSetParser()
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(
            options.output_dir, binIds, options.marker_file)

        # compare markers found in each bin

        merger = Merger()
        outputFile = merger.run(binFiles, options.output_dir, "merger.table.txt", binIdToModels, binIdToBinMarkerSets,
                                options.delta_comp, options.delta_cont, options.merged_comp, options.merged_cont)

        self.logger.info('Merger information written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def outliers(self, options):
        """Outlier command"""

        self.logger.info('[CheckM - outlier] Identifying outliers in bins.')

        checkBinInputExists(options.bin_input, False)
        checkFileExists(options.tetra_profile)
        makeSurePathExists(os.path.dirname(options.output_file))

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        binTools = BinTools()
        binTools.identifyOutliers(options.results_dir,
                                  binFiles,
                                  options.tetra_profile,
                                  options.distributions,
                                  options.report_type,
                                  options.output_file)

        self.logger.info('Outlier information written to: ' +
                         options.output_file)

        self.timeKeeper.printTimeStamp()

    def modify(self, options):
        """Modify command"""

        self.logger.info('[CheckM - modify] Modifying sequences in bin.')

        makeSurePathExists(os.path.dirname(options.output_file))

        if not (options.add or options.remove or options.outlier_file):
            self.logger.error('No modification to bin requested.\n')
            sys.exit(1)

        if (options.add or options.remove) and options.outlier_file:
            self.logger.error(
                "The 'outlier_file' option cannot be specified with 'add' or 'remove'.\n")
            sys.exit(1)

        binTools = BinTools()

        if options.add or options.remove:
            binTools.modify(options.bin_file, options.seq_file,
                            options.add, options.remove, options.output_file)
        elif options.outlier_file:
            binTools.removeOutliers(
                options.bin_file, options.outlier_file, options.output_file)

        self.logger.info('Modified bin written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def unique(self, options):
        """Unique command"""

        self.logger.info(
            '[CheckM - unique] Ensuring no sequences are assigned to multiple bins.')

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        binTools = BinTools()
        binTools.unique(binFiles)

        self.timeKeeper.printTimeStamp()

    def ssuFinder(self, options):
        """SSU finder command"""

        self.logger.info(
            '[CheckM - ssu_finder] Identifying SSU (16S/18S) rRNAs in sequences.')

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        checkFileExists(options.seq_file)
        makeSurePathExists(options.output_dir)

        ssuFinder = SSU_Finder(options.threads)
        ssuFinder.run(options.seq_file, binFiles, options.output_dir,
                      options.evalue, options.concatenate)

        self.timeKeeper.printTimeStamp()

    def test(self, options):
        """Quick test of CheckM"""

        self.logger.info(
            '[CheckM - Test] Processing E.coli K12-W3310 to verify operation of CheckM.')

        verifyEcoli = VerifyEcoli()
        verifyEcoli.run(self, options.output_dir)

        self.timeKeeper.printTimeStamp()

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""

        try:
            if options.file == "stdout":
                options.file = ''
        except:
            pass

        if(options.subparser_name == "data"):
            self.updateCheckM_DB(options)
        elif(options.subparser_name == 'tree'):
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
        elif(options.subparser_name == 'lineage_wf'):
            options.marker_file = os.path.join(
                options.output_dir, 'lineage.ms')
            options.tree_dir = options.output_dir
            options.analyze_dir = options.output_dir
            options.out_format = 1
            options.bAlignTopHit = False
            options.exclude_markers = None
            options.coverage_file = None

            self.tree(options)
            self.lineageSet(options)
            self.analyze(options)
            self.qa(options)
        elif(options.subparser_name == 'taxonomy_wf'):
            options.marker_file = os.path.join(
                options.output_dir, options.taxon + '.ms')
            options.analyze_dir = options.output_dir
            options.out_format = 1
            options.bAlignTopHit = False
            options.exclude_markers = None

            self.taxonSet(options)
            self.analyze(options)
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
        elif(options.subparser_name == 'len_hist'):
            self.lengthHistogram(options)
        elif(options.subparser_name == 'marker_plot'):
            self.markerPlot(options)
        elif(options.subparser_name == 'gc_bias_plot'):
            self.gcBiasPlot(options)
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
        elif(options.subparser_name == 'modify'):
            self.modify(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'ssu_finder'):
            self.ssuFinder(options)
        elif(options.subparser_name == 'test'):
            options.bCalledGenes = False
            options.pplacer_threads = 1
            self.test(options)
        else:
            self.logger.error('Unknown CheckM command: ' +
                              options.subparser_name + '\n')
            sys.exit(1)

        return 0
