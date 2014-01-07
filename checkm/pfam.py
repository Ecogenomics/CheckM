###############################################################################
#
# pfam.py - methods for identifying PFAM genes from HMMER models.
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
import operator

from common import checkFileExists

class PFAM(object):
    def __init__(self):
        self.pfamClanFile = os.path.join(os.path.dirname(sys.argv[0]), '..', 'data', 'pfam', 'Pfam-A.hmm.dat')
        self.idToAcc = {}   # map ids to pfam accessions
        self.clan = {}      # map pfam accessions to clans
        self.nested = {}    # map pfam accessions to nested pfam accessions
    
    def __readClansAndNesting(self):     
        checkFileExists(self.pfamClanFile)
     
        idNested = {}
        for line in open(self.pfamClanFile):
            if '#=GF ID' in line:
                ID = line.split()[2].strip()
            elif '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
                pfamAcc = pfamAcc.replace('PF', 'pfam')
                pfamAcc = pfamAcc[0:pfamAcc.rfind('.')]
                self.idToAcc[ID] = pfamAcc
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                self.clan[pfamAcc] = clanId
            elif '#=GF NE' in line:
                nestedId = line.split()[2].strip()
    
                if nestedId not in idNested:
                    idNested[nestedId] = []
                idNested[nestedId].append(ID)
    
                if ID not in idNested:
                    idNested[ID] = []
                idNested[ID].append(nestedId)
    
        # set nested structure to use pfam accessions instead of IDs
        for ID, nested in idNested.iteritems():
            pfamAcc = self.idToAcc[ID]
            self.nested[pfamAcc] = set([self.idToAcc[x] for x in nested])
    
    def pfamIdToClanId(self):
        checkFileExists(self.pfamClanFile)
        
        d = {}
        for line in open(self.pfamClanFile):
            if '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
                pfamAcc = pfamAcc.replace('PF', 'pfam')
                pfamAcc = pfamAcc[0:pfamAcc.rfind('.')]
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                d[pfamAcc] = clanId
    
        return d
    
    def filterHitsFromSameClan(self, pfamHits, pfamMarkers):
        # check if clan and nesting need to be computer
        if len(self.clan) == 0:
            self.__readClansAndNesting()
    
        # for each gene, take only the best hit for each PFAM clan
        hitsToPfamMarkers = {}
        for geneId, hits in pfamHits.iteritems():
            # sort in ascending order of evalue
            hits.sort(key=operator.itemgetter(1))
    
            filtered = set()
            for i in xrange(0, len(hits)):
                if i in filtered:
                    continue
    
                pfamIdI = hits[i][0]
                clanI = self.clan.get(pfamIdI, None)
                startI = hits[i][2]
                endI = hits[i][3]
    
                for j in xrange(i+1, len(hits)):
                    if j in filtered:
                        continue
    
                    pfamIdJ = hits[j][0]
                    clanJ = self.clan.get(pfamIdJ, None)
                    startJ = hits[j][2]
                    endJ = hits[j][3]
    
                    # check if hits are from the same clan
                    if pfamIdI != None and pfamIdJ != None and clanI == clanJ:
                        # check if hits overlap
                        if (startI <= startJ and endI > startJ) or (startJ <= startI and endJ > startI):
                            # check if pfams are nested
                            if not (pfamIdI in self.nested and pfamIdJ in self.nested[pfamIdI]):
                                # hits should be filtered as it is from the same clan, overlaps, and is not
                                # nested with a pfam hit with a lower e-value
                                filtered.add(j)
    
            # tabulate unfiltered hits
            for i in xrange(0, len(hits)):
                if i in filtered:
                    continue
    
                pfamId = hits[i][0]
                if pfamId in pfamMarkers:
                    s = hitsToPfamMarkers.get(pfamId, set())
                    s.add(geneId)
                    hitsToPfamMarkers[pfamId] = s
    
        return hitsToPfamMarkers
