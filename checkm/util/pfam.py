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

from collections import defaultdict

from checkm.common import checkFileExists


def stripVersion(accession):
    """Remove version suffix (e.g., '.14' in 'PF00162.14') from a Pfam accession, if present."""
    if '.' in accession:
        return accession[0:accession.rfind('.')]

    return accession


class PFAM(object):
    def __init__(self, pfamClanFile):
        self.pfamClanFile = pfamClanFile
        self.idToAcc = {}  # map ids to pfam accessions
        self.clan = {}  # map pfam accessions to clans
        self.nested = {}  # map pfam accessions to nested pfam accessions

    def __readClansAndNesting(self):
        checkFileExists(self.pfamClanFile)

        idNested = defaultdict(list)
        for line in open(self.pfamClanFile):
            if '#=GF ID' in line:
                ID = line.split()[2].strip()
            elif '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
                pfamAcc = stripVersion(pfamAcc)
                self.idToAcc[ID] = pfamAcc
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                self.clan[pfamAcc] = clanId
            elif '#=GF NE' in line:
                nestedId = line.split()[2].strip()
                idNested[nestedId].append(ID)
                idNested[ID].append(nestedId)

        # set nested structure to use pfam accessions instead of IDs
        for ID, nested in idNested.items():
            pfamAcc = self.idToAcc[ID]
            self.nested[pfamAcc] = set([self.idToAcc[x] for x in nested])

    def pfamIdToClanId(self):
        """Determine clan of each pfam."""
        checkFileExists(self.pfamClanFile)

        d = {}
        for line in open(self.pfamClanFile):
            if '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                d[pfamAcc] = clanId

        return d

    def genesInClan(self):
        """Determine all genes within each clan."""
        checkFileExists(self.pfamClanFile)

        d = defaultdict(set)
        for line in open(self.pfamClanFile):
            if '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                d[clanId].update([pfamAcc])

        return d

    def filterHitsFromSameClan(self, markerHits):
        """Filter hits to ORF from same clan."""

        # check if clan and nesting need to be computed
        if len(self.clan) == 0:
            self.__readClansAndNesting()

        # determine all PFAM hits to each ORF and setup initial set of filtered markers
        # (a plain dict is used so looking up an absent marker never inserts an empty hit list)
        filteredMarkers = {}
        hitsToORFs = defaultdict(list)
        for markerId, hits in markerHits.items():
            if markerId.startswith('PF'):
                for hit in hits:
                    hitsToORFs[hit.target_name].append(hit)
            else:
                # retain all non-PFAM markers
                filteredMarkers[markerId] = hits

        # for each gene, take only the best hit for each PFAM clan
        for target_name, hits in hitsToORFs.items():
            # sort in ascending order of full sequence e-value followed by domain (independent) e-value
            hits.sort(key=lambda x: (x.full_e_value, x.i_evalue))

            filtered = set()
            for i in range(0, len(hits)):
                if i in filtered:
                    continue

                pfamIdI = stripVersion(hits[i].query_accession)
                clanI = self.clan.get(pfamIdI, None)
                startI = hits[i].ali_from
                endI = hits[i].ali_to

                for j in range(i + 1, len(hits)):
                    if j in filtered:
                        continue

                    pfamIdJ = stripVersion(hits[j].query_accession)
                    clanJ = self.clan.get(pfamIdJ, None)
                    startJ = hits[j].ali_from
                    endJ = hits[j].ali_to

                    # check if hits are from the same clan (Pfams without a clan are unrelated)
                    if clanI is not None and clanI == clanJ:
                        # check if hits overlap
                        if (startI <= startJ and endI > startJ) or (startJ <= startI and endJ > startI):
                            # check if pfams are nested
                            if not (pfamIdI in self.nested and pfamIdJ in self.nested[pfamIdI]):
                                # hit should be filtered as it is from the same clan, overlaps, and is not
                                # nested with a pfam hit with a lower e-value
                                filtered.add(j)

            # tabulate unfiltered hits
            for i in range(0, len(hits)):
                if i in filtered:
                    continue

                filteredMarkers.setdefault(hits[i].query_accession, []).append(hits[i])

        return filteredMarkers

    def genesInSameClan(self, genes):
        """Get all genes from the PFAM clans spanned by the input gene set."""

        # get a list of clans spanned by the input gene set
        pfamIdToClanId = self.pfamIdToClanId()

        clans = set()
        for gene in genes:
            clanId = pfamIdToClanId.get(gene, None)
            if clanId != None:
                clans.add(clanId)

        # get a list of all other genes from these clans
        genesInClan = self.genesInClan()

        allGenesInClans = set()
        for clan in clans:
            allGenesInClans.update(genesInClan[clan])

        return allGenesInClans - genes
