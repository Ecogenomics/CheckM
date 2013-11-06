#!/usr/bin/env python

###############################################################################
#
# runDistDeltaGC.py - calculate delta GC distribution of reference genomes
#                     in parallel
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

if __name__ == '__main__':
    numProcessors = 24

    startGC = 0.0
    endGC = 1.0
    gcStep = 0.01

    curGC = startGC
    cmdList = open('cmdList.txt', 'w')
    while curGC < endGC:
        cmdList.write('python distributionDeltaGC.py --len_step_size 200 --max_len 5000 --gc_start ' + str(curGC) + ' --gc_end ' + str(curGC + gcStep) + ' /srv/db/img/4.1/metadata/img_metadata_4_1.dhp.tsv  /srv/db/img/4.1/genomes/' + '\n')
        curGC += gcStep

    cmdList.close()

    os.system('cat cmdList.txt | parallel --max-procs ' + str(numProcessors))
