#!/usr/bin/env python
###############################################################################
#                                                                             #
#    metachecka2000                                                           #
#                                                                             #
#    Entry point. See metachecka2000/metachecka2000.py for internals          #
#                                                                             #
#    Copyright (C) Michael Imelfort                                     #
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

import argparse
import sys
import re
from metachecka2000 import metachecka2000 

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def printHelp():
    print '''\
    
             ...::: metachecka2000 :::...
                   
    Betta check your meta before you wreck your meta!

    metachecka2000 build   -> Create the underlying data needed for testing
    metachecka2000 qa      -> Test for bin contamination and completeness
    metachecka2000 all     -> Build and test in one step
    
    USE: metachecka2000 OPTION -h to see detailed options
    '''

if __name__ == '__main__':

    #-------------------------------------------------
    # intialise the options parser
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    ##################################################
    # Typical workflow
    ##################################################
    
    #-------------------------------------------------
    # parse raw data and save
    file_parser = subparsers.add_parser('parse',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='parse raw data and save to disk',
                                        description='Parse raw data and save to disk')
    file_parser.add_argument('dbname', help="name of the database being created")
    file_parser.add_argument('reference', help="fasta file containing bam reference sequences")
    file_parser.add_argument('bamfiles', nargs='+', help="bam files to parse")
    file_parser.add_argument('-d', '--dump', action="store_true", default=False, help="dump the contents of all tables to screen")
    file_parser.add_argument('-f', '--force', action="store_true", default=False, help="overwrite existing DB file without prompting")
        
        
    #-------------------------------------------------
    # get and check options
    args = None
    if(len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
        printHelp()
        sys.exit(0)
    else:
        args = parser.parse_args()

    #-------------------------------------------------
    # do what we came here to do
    try:
        MC2K_parser = groopm.MC2KOptionsParser()
        if(False):
            #import pstats
            #p = pstats.Stats('prof')
            #p.sort_stats('cumulative').print_stats(10)
            #p.sort_stats('time').print_stats(10)
            import cProfile
            cProfile.run('MC2K_parser.parseOptions(args)', 'prof')
        else:        
            MC2K_parser.parseOptions(args)
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
        
