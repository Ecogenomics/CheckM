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

import resultsParser
import markerGeneFinder
import database as chmdb

def connectToDatabase(database_name):
    """ Return a database object based on name """
    if database_name is None:
        database_name = os.getenv("CHECKM_DB")
        if database_name is None:
            raise RuntimeError('Cannot connect to DB')

    return chmdb.MarkerDB(db=database_name)


class OptionsParser():
    def __init__(self): pass

    def build(self, options, db=None):
        """Build command"""
        if not options.quiet:
            print "Building data prior to checking."
                
        mgf = markerGeneFinder.MarkerGeneFinder(threads=options.threads)
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
        mgf.find(target_files, options.out_folder, options.hmm, quiet=options.quiet)

    def qa(self, options):
        """QA Command"""
        if not options.quiet:
            print "Analysing bins."
                
        RP = resultsParser.HmmerResultsParser()
        RP.analyseResults(options.out_folder,
                          options.hmm,
                          eCO=options.e_value,
                          lengthCO=options.length,
                          quiet=options.quiet,
                          outFile=options.file
                          )
        RP.printSummary(outputFormat=options.out_format)

    def align(self, options):
        """Align Command"""
        if not options.quiet:
            print "Constructing alignments."
                
        if hasattr(options, 'separate'):
            HA = resultsParser.HMMAligner(options.separate, options.consensus, options.out_format)
        else:
            HA = resultsParser.HMMAligner()
        bh = False
        if hasattr(options, 'best_alignment'):
            bh=True
        HA.makeAlignments(options.out_folder,
                          options.hmm,
                          eCO=options.e_value,
                          lengthCO=options.length,
                          quiet=options.quiet,
                          bestHit=bh
                          )

    def makeDB(self, options):
        DB = chmdb.MarkerDB()
        DB.makeDB(options)

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""
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

        if(options.subparser_name == 'build'):
            self.build(options)

        elif(options.subparser_name == 'qa'):
            self.qa(options)

        elif(options.subparser_name == 'all'):
            self.build(options)
            self.qa(options)
            self.align(options)

        elif(options.subparser_name == 'align'):
            self.align(options)

        if database_name is not None:
            os.remove(tmp)

        return 0
