#!/usr/bin/env python

###############################################################################
#
# runGenomeSigDist.py - calculate distribution of tetranucleotide
#                             distances for all reference genomes
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
import argparse

class RunGenomeSigDist(object):
    def __init__(self):
        pass

    def run(self, metadataFile, genomeDir, numWindows, numThreads):
        # read metadata file
        print 'Determining finished prokaryotic reference genomes.'
        genomeIds = []
        bHeader = True
        for line in open(metadataFile):
            lineSplit = [x.strip() for x in line.split('\t')]

            if bHeader:
                bHeader = False
                continue

            domain = lineSplit[1]
            status = lineSplit[2]
            if status == 'Finished' and (domain == 'Bacteria' or domain == 'Archaea'):
                genomeId = lineSplit[0]

                if os.path.exists(os.path.join(genomeDir, genomeId, genomeId + '.fna')):
                    genomeIds.append(genomeId)

        print 'Finished genomes to process: ' + str(len(genomeIds))

        # calculate difference in genomic signatures for each genome
        fout = open('cmdList.txt', 'w')
        for genomeId in genomeIds:
            genomeFile = os.path.join(genomeDir, genomeId, genomeId + '.fna')
            fout.write('./GenomeSigDist/genome-sig-dist -n ' + str(numWindows) + ' -g ' + genomeFile + ' -o ./deltaTD/' + genomeId + '.txt' + '\n')
        fout.close()

        os.system('cat cmdList.txt | parallel --max-procs ' + str(numThreads))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate GC distribution over reference genomes.')
    parser.add_argument('metadata_file', help='IMG metadata file.')
    parser.add_argument('genome_dir', help='IMG genome directory.')
    parser.add_argument('--num_windows', help='number of windows to sample', type = int, default = 10000)
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)

    args = parser.parse_args()

    runGenomeSigDist = RunGenomeSigDist()
    runGenomeSigDist.run(args.metadata_file, args.genome_dir, args.num_windows, args.threads)
