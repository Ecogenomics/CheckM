#!/usr/bin/env python
###############################################################################
#                                                                             #
#    metachecka2000                                                           #
#                                                                             #
#    Wraps coarse workflows                                                   #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#                888b     d888  .d8888b.   .d8888b.  888    d8P               #  
#                8888b   d8888 d88P  Y88b d88P  Y88b 888   d8P                #  
#                88888b.d88888 888    888        888 888  d8P                 #
#                888Y88888P888 888             .d88P 888d88K                  #
#                888 Y888P 888 888         .od888P"  8888888b                 #
#                888  Y8P  888 888    888 d88P"      888  Y88b                #
#                888   "   888 Y88b  d88P 888"       888   Y88b               #
#                888       888  "Y8888P"  888888888  888    Y88b              #
#                                                                             #
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

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
import os 
import argparse

# MetaChecka2000 imports
import resultsParser
import dataConstructor

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Mc2kOptionsParser():
    def __init__(self): pass
    
    def Mc2kBuild(self, options):
        """Build command"""
        DC = dataConstructor.Mc2kHmmerDataConstructor(threads=options.threads)
        target_files = []
        if options.bin_folder is not None:
            all_files = os.listdir(options.bin_folder)
            for j in all_files:
                if j.endswith(options.extension):
                    target_files.append(os.path.join(options.bin_folder,j))
        
        if options.infiles:
            target_files.extend(options.infiles)
        
        if not target_files:
            raise "No input files!"
        DC.buildData(target_files,
                     options.out_folder,
                     options.hmm,
                     options.closed,
                     options.prefix,
                     verbose=options.verbose
                     )
    
    def Mc2kQa(self, options):
        """QA Command"""
        RP = resultsParser.Mc2kHmmerResultsParser(prefix=options.prefix)
        RP.analyseResults(options.out_folder,
                          options.hmm,
                          eCO=options.e_value,
                          lengthCO=options.length,
                          verbose=options.verbose,
                          outFile=options.file
                          )

    def Mc2kAlign(self, options):
        """Align Command"""
        if hasattr(options, 'separate'):
            HA = resultsParser.HMMAligner(options.prefix,
                      options.separate,
                      options.consensus,
                      options.out_format
                      )
        else:
            HA = resultsParser.HMMAligner(options.prefix)
        bh = False
        if hasattr(options, 'best_alignment'):
            bh=True
        HA.makeAlignments(options.out_folder,
                          options.hmm,
                          eCO=options.e_value,
                          lengthCO=options.length,                              
                          verbose=options.verbose,
                          prefix=options.prefix,
                          bestHit=bh
                          )

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""
        try:
            if options.file == "STDOUT":
                options.file = ''
        except:
            pass

        if(options.subparser_name == 'build'):
            # build prodigal and hmm result
            if options.verbose:
                print "Building data prior to checking..."
            self.Mc2kBuild(options)
                            
        elif(options.subparser_name == 'qa'):
            # do qa analysis
            if options.verbose:
                print "Analysing bins..."
            self.Mc2kQa(options)

        elif(options.subparser_name == 'all'):
            # all in one
            print "Building data prior to checking..."
            self.Mc2kBuild(options)
            print "Analysing bins..."
            self.Mc2kQa(options)
            print "Constructing alignments..."
            self.Mc2kAlign(options)
            
        elif(options.subparser_name == 'align'):
            
            # make alignments
            if options.verbose:
                print "Constructing alignments..."
            self.Mc2kAlign(options)

        return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################
