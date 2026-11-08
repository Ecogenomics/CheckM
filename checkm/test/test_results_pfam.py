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
import shutil
import tempfile
import unittest

from checkm.util.pfam import PFAM


class Hit(object):
    """Minimal stand-in for checkm.hmmer.HmmerHitDOM."""

    def __init__(self, target_name, query_accession, full_e_value, ali_from, ali_to, i_evalue=None):
        self.target_name = target_name
        self.query_accession = query_accession
        self.full_e_value = full_e_value
        self.i_evalue = full_e_value if i_evalue is None else i_evalue
        self.ali_from = ali_from
        self.ali_to = ali_to


PFAM_DAT = """# STOCKHOLM 1.0
#=GF ID   FamA
#=GF AC   PF00001.21
#=GF CL   CL0001
//
#=GF ID   FamB
#=GF AC   PF00002.9
#=GF CL   CL0001
#=GF NE   FamC
//
#=GF ID   FamC
#=GF AC   PF00003.5
#=GF CL   CL0001
#=GF NE   FamB
//
#=GF ID   FamD
#=GF AC   PF00004.1
//
#=GF ID   FamE
#=GF AC   PF00005.2
//
#=GF ID   FamF
#=GF AC   PF00006.3
#=GF CL   CL0002
//
"""


class VerifyPfamClanFilter(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        self.datFile = os.path.join(self.tmpDir, 'Pfam-A.hmm.dat')
        with open(self.datFile, 'w') as f:
            f.write(PFAM_DAT)
        self.pfam = PFAM(self.datFile)

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def _filter(self, hits):
        markerHits = {}
        for h in hits:
            markerHits.setdefault(h.query_accession, []).append(h)
        return self.pfam.filterHitsFromSameClan(markerHits)

    def testSameClanOverlapFiltered(self):
        """Overlapping hits to one ORF from the same clan keep only the best hit."""
        best = Hit('c_1', 'PF00001.21', 1e-50, 10, 100)
        worse = Hit('c_1', 'PF00002.9', 1e-20, 50, 150)
        filtered = self._filter([worse, best])
        self.assertEqual(filtered, {'PF00001.21': [best]})

    def testSameClanNoOverlapRetained(self):
        a = Hit('c_1', 'PF00001.21', 1e-50, 10, 100)
        b = Hit('c_1', 'PF00002.9', 1e-20, 100, 150)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00001.21': [a], 'PF00002.9': [b]})

    def testSameClanDifferentOrfRetained(self):
        a = Hit('c_1', 'PF00001.21', 1e-50, 10, 100)
        b = Hit('c_2', 'PF00002.9', 1e-20, 10, 100)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00001.21': [a], 'PF00002.9': [b]})

    def testNestedRetained(self):
        a = Hit('c_1', 'PF00002.9', 1e-50, 10, 100)
        b = Hit('c_1', 'PF00003.5', 1e-20, 20, 80)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00002.9': [a], 'PF00003.5': [b]})

    def testClanlessOverlapRetained(self):
        """Two overlapping Pfams without a clan are unrelated and must both be kept."""
        a = Hit('c_1', 'PF00004.1', 1e-50, 10, 100)
        b = Hit('c_1', 'PF00005.2', 1e-20, 20, 120)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00004.1': [a], 'PF00005.2': [b]})

    def testClanlessAndClanOverlapRetained(self):
        a = Hit('c_1', 'PF00004.1', 1e-50, 10, 100)
        b = Hit('c_1', 'PF00001.21', 1e-20, 20, 120)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00004.1': [a], 'PF00001.21': [b]})

    def testDifferentClanOverlapRetained(self):
        a = Hit('c_1', 'PF00006.3', 1e-50, 10, 100)
        b = Hit('c_1', 'PF00001.21', 1e-20, 20, 120)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00006.3': [a], 'PF00001.21': [b]})

    def testAccessionWithoutVersion(self):
        """Accessions without a version suffix must not lose their last character."""
        best = Hit('c_1', 'PF00001', 1e-50, 10, 100)
        worse = Hit('c_1', 'PF00002', 1e-20, 50, 150)
        filtered = self._filter([worse, best])
        self.assertEqual(filtered, {'PF00001': [best]})

        # 'PF00010' would become 'PF0001' if the last character was stripped,
        # which must not be confused with a clan member
        a = Hit('c_2', 'PF00010', 1e-50, 10, 100)
        b = Hit('c_2', 'PF00001.21', 1e-20, 20, 120)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'PF00010': [a], 'PF00001.21': [b]})

    def testNonPfamRetained(self):
        a = Hit('c_1', 'TIGR00001', 1e-50, 10, 100)
        b = Hit('c_1', 'TIGR00002', 1e-20, 10, 100)
        filtered = self._filter([a, b])
        self.assertEqual(filtered, {'TIGR00001': [a], 'TIGR00002': [b]})

    def testReturnsPlainDict(self):
        """Looking up an absent marker must not insert an empty hit list."""
        filtered = self._filter([Hit('c_1', 'PF00001.21', 1e-50, 10, 100)])
        self.assertNotIn('PF00005.2', filtered)
        self.assertRaises(KeyError, lambda: filtered['PF00005.2'])
        self.assertNotIn('PF00005.2', filtered)


if __name__ == "__main__":
    unittest.main()
