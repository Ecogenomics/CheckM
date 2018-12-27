#!/usr/bin/env python

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

"""
Get all PFAM and TIGRFAM HMMs required for marker genes.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import argparse

from lib.img import IMG
from lib.pfam import PFAM
from lib.markerSet import MarkerSet

class GetHMMs(object):
    def __init__(self):
        pass

    def run(self, minGenomes, minMarkerSets):
        img = IMG()
        pfam = PFAM()

        # get list of all marker genes
        markerset = MarkerSet()
        pfamIds, tigrIds = markerset.getCalculatedMarkerGenes()

        print('TIGR marker genes: ' + str(len(tigrIds)))
        print('PFAM marker genes: ' + str(len(pfamIds)))

        # get all PFAM HMMs that are in the same clan
        # as any of the marker genes
        pfamIdToClanId = pfam.pfamIdToClanId()
        clans = set()
        for pfamId in pfamIds:
            if pfamId.replace('PF', 'pfam') in pfamIdToClanId:
                clans.add(pfamIdToClanId[pfamId])

        for pfamId, clanId in pfamIdToClanId.items():
            if clanId in clans:
                pfamIds.add(pfamId)

        print('  PFAM HMMs require to cover marker gene clans: ' + str(len(pfamIds)))

        # get name of each PFAM HMM
        fout = open('./hmm/pfam.keyfile.txt', 'w')
        pfamNames = []
        for line in open(img.pfamHMMs):
            if 'NAME' in line:
                name = line[line.find(' '):].strip()
            elif 'ACC' in line:
                acc = line[line.find(' '):line.rfind('.')].strip()
                if acc.replace('PF', 'pfam') in pfamIds:
                    pfamNames.append(name)
                    fout.write(name + '\n')
        fout.close()

        print('PFAM names: ' + str(len(pfamNames)))

        # extract each PFAM HMM
        os.system('hmmfetch -f ' + img.pfamHMMs + ' ./hmm/pfam.keyfile.txt > ./hmm/pfam_markers.hmm')

        # get name of each PFAM HMM
        fout = open('./hmm/tigr.keyfile.txt', 'w')
        for tigrId in tigrIds:
            fout.write(tigrId + '\n')
        fout.close()

        # extract each PFAM HMM
        os.system('hmmfetch -f ' + img.tigrHMMs + ' ./hmm/tigr.keyfile.txt > ./hmm/tigr_markers.hmm')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get all PFAM HMMs.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-m', '--min_genomes', help='Minimum genomes required to include in analysis', type=int, default = 20)
    parser.add_argument('-x', '--min_markers', help='Minimum marker sets required to include in analysis', type=int, default = 20)

    args = parser.parse_args()

    getHMMs = GetHMMs()
    getHMMs.run(args.min_genomes, args.min_markers)
