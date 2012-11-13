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
__version__ = "0.0.1"
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

class MC2KOptionsParser():
    def __init__(self): pass
    
    def parseOptions(self, options ):

        if(options.subparser_name == 'build'):
            # build prodigal and hmm result
            if options.verbose:
                print "Building data prior to checking..."
            DC = dataConstructor.MC2KHmmerDataConstructor(threads=options.threads)
            DC.buildData(options.bin_folder,
                         options.out_folder,
                         options.hmm,
                         options.extension,
                         options.closed,
                         options.prefix,
                         verbose=options.verbose
                         )
                            
        elif(options.subparser_name == 'qa'):
            # do qa analysis
            if options.verbose:
                print "Checking bins..."
            RP = resultsParser.MC2KHmmerResultsParser(prefix=options.prefix)
            RP.analyseResults(options.out_folder,
                              options.hmm,
                              eCO=options.e_value,
                              lengthCO=options.length,
                              verbose=options.verbose,
                              outFile=options.file
                              )

        elif(options.subparser_name == 'all'):
            # all in one
            if options.verbose:
                print "Running complete MC2K pipeline..."
            if options.verbose:
                print "Building data prior to checking..."
            DC = dataConstructor.MC2KHmmerDataConstructor(threads=options.threads)
            DC.buildData(options.bin_folder,
                         options.out_folder,
                         options.hmm,
                         options.extension,
                         options.closed,
                         options.prefix,
                         verbose=options.verbose
                         )
            if options.verbose:
                print "Checking bins..."
            RP = resultsParser.MC2KHmmerResultsParser(prefix=options.prefix)
            RP.analyseResults(options.out_folder,
                              options.hmm,
                              eCO=options.e_value,
                              lengthCO=options.length,
                              verbose=options.verbose,
                              outFile=options.file
                              )
        return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################
