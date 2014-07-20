#!/usr/bin/env python

###############################################################################
#
# tigr2pfam.py - find TIGRFAMs that map to the same gene as a PFAM
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

__prog_desc__ ='find TIGRFAMs that map to the same gene as a PFAM'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import os
import argparse

from checkm.lib.img import IMG

class Tigr2Pfam(object):
    def __init__(self):
        pass

    def run(self, metadataFile, percentThreshold):
        img = IMG()

        metadata = img.genomeMetadataFromFile(metadataFile)

        matches = {}
        pfamCount = {}
        tigrCount = {}
        for genomeCounter, genomeId in enumerate(metadata):
            statusStr = '  Finished processing %d of %d (%.2f%%) genomes.' % (genomeCounter+1, len(metadata), float(genomeCounter+1)*100/len(metadata))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            if metadata[genomeId]['status'] == 'Finished':
                pfamFile = img.genomeDir + genomeId + '/' + genomeId + img.pfamExtension

                if not os.path.exists(pfamFile):
                    continue

                # get PFAM hits
                geneIdToPfams = {}
                bHeader = True
                for line in open(pfamFile):
                    if bHeader:
                        bHeader = False
                        continue

                    lineSplit = line.split('\t')
                    if lineSplit[0] in geneIdToPfams:
                        geneIdToPfams[lineSplit[0]].add(lineSplit[8])
                    else:
                        geneIdToPfams[lineSplit[0]] = set([lineSplit[8]])

                    if lineSplit[8] in pfamCount:
                        pfamCount[lineSplit[8]].add(genomeId)
                    else:
                        pfamCount[lineSplit[8]] = set([genomeId])

                # get TIGRFAM hits
                geneIdToTigr = {}
                bHeader = True
                for line in open(img.genomeDir + genomeId + '/' + genomeId + img.tigrExtension):
                    if bHeader:
                        bHeader = False
                        continue

                    lineSplit = line.split('\t')
                    if lineSplit[0] in geneIdToTigr:
                        geneIdToTigr[lineSplit[0]].add(lineSplit[6])
                    else:
                        geneIdToTigr[lineSplit[0]] = set([lineSplit[6]])

                    if lineSplit[6] in tigrCount:
                        tigrCount[lineSplit[6]].add(genomeId)
                    else:
                        tigrCount[lineSplit[6]] = set([genomeId])

                # keep track of TIGRFAMs matching the same gene as a PFAM
                geneIds = set(geneIdToPfams.keys()).union(set(geneIdToTigr.keys()))
                for geneId in geneIds:
                    pfams = geneIdToPfams.get(geneId, None)
                    tigrs = geneIdToTigr.get(geneId, None)

                    if pfams == None or tigrs == None:
                        continue

                    for pfamId in pfams:
                        for tigrId in tigrs:
                            key = pfamId + '-' + tigrId
                            if key in matches:
                                matches[key].add(genomeId)
                            else:
                                matches[key] = set([genomeId])

        sys.stdout.write('\n')

        # find TIGRFAMs that generally hit the same gene as a PFAM
        fout = open('../data/pfam/tigrfam2pfam.tsv', 'w')
        for key, genomeSet in matches.iteritems():
            pfam, tigr = key.split('-')

            # deem a TIGRFAM HMM redundant if it is almost always hits that
            # same ORF as a PFAM HMM
            if float(len(genomeSet)) / len(tigrCount[tigr]) >= percentThreshold:
                fout.write(pfam + '\t' + tigr + '\n')
        fout.close()

if __name__ == '__main__':
    print 'Tigr2Pfam v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('metadata_file', help='IMG metadata file')
    parser.add_argument('-p', '--percent', help='percent threshold for matching a TIGRFAM to a PFAM', type = float, default = 0.90)

    args = parser.parse_args()

    tigr2Pfam = Tigr2Pfam()
    tigr2Pfam.run(args.metadata_file, args.percent)
