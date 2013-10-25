###############################################################################
#
# hmmer.py - runs HMMER and provides functions for parsing output  
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

import defaultValues

from common import makeSurePathExists

class ProdigalError(BaseException): pass

class ProdigalRunner():
    """Wrapper for running HMMER3"""
    def __init__(self):
        # make sure prodigal is installed
        self.checkForProdigal()
        
    def run(self, query, outputDir, translationTable=11):
        makeSurePathExists(outputDir)
        
        aaFile = os.path.join(outputDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_AA__)
        ntFile = os.path.join(outputDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_NT__)
        gffFile = os.path.join(outputDir, defaultValues.__CHECKM_DEFAULT_PRODIGAL_GFF__)
        
        cmd = ('prodigal -q -c -m -f gff -t %s -a %s -d %s -i %s > %s' % (str(translationTable), aaFile, ntFile, query, gffFile))
        os.system(cmd)
        
        return aaFile

    def checkForProdigal(self):
        """Check to see if Prodigal is on the system before we try to run it.
    
        We assume that a successful prodigal -h returns 0 and anything
        else returns something non-zero
        """
        # redirect stdout so we don't get mess!
        try:
            exit_status = os.system('prodigal -h 2> /dev/null')
        except:
            print "Unexpected error!", sys.exc_info()[0]
            raise
    
        if exit_status != 0:
            raise ProdigalError("Error attempting to run prodigal, is it in your path?")

class ProdigalParser():
    """Parses tabular output"""
    def __init__(self, fileHandle):
        """Give this guy an open file!"""
        self.handle = fileHandle

    def next(self):
        """Get the next result in the file"""
        pass
