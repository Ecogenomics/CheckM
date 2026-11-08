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

import unittest

from checkm.resultsParser import ResultsManager


class Hit(object):
    """Minimal stand-in for checkm.hmmer.HmmerHitDOM."""

    def __init__(self, target_name, query_accession,
                 hmm_from, hmm_to, ali_from, ali_to, env_from, env_to,
                 target_length=200, full_e_value=1e-30):
        self.target_name = target_name
        self.query_accession = query_accession
        self.target_length = target_length
        self.full_e_value = full_e_value
        self.hmm_from = hmm_from
        self.hmm_to = hmm_to
        self.ali_from = ali_from
        self.ali_to = ali_to
        self.env_from = env_from
        self.env_to = env_to


class VerifyAdjacentMarkerGenes(unittest.TestCase):
    def _manager(self, markerHits):
        rm = ResultsManager('bin1', {})
        rm.markerHits = markerHits
        return rm

    def testMergeCoordinates(self):
        """Merged hits span from the smallest start to the largest end coordinate."""
        a = Hit('contig_1_5', 'PF00001.1', 1, 120, 3, 118, 1, 120, target_length=130)
        b = Hit('contig_1_6', 'PF00001.1', 110, 300, 2, 190, 1, 195, target_length=210)
        rm = self._manager({'PF00001.1': [b, a]})
        rm.identifyAdjacentMarkerGenes()

        hits = rm.markerHits['PF00001.1']
        self.assertEqual(len(hits), 1)
        h = hits[0]
        self.assertEqual(h.target_name, 'contig_1_5&&contig_1_6')
        self.assertEqual(h.target_length, 340)
        self.assertEqual((h.hmm_from, h.hmm_to), (1, 300))
        self.assertEqual((h.ali_from, h.ali_to), (2, 190))
        self.assertEqual((h.env_from, h.env_to), (1, 195))

    def testNonAdjacentNotMerged(self):
        a = Hit('contig_1_5', 'PF00001.1', 1, 120, 3, 118, 1, 120)
        b = Hit('contig_1_7', 'PF00001.1', 110, 300, 2, 190, 1, 195)
        c = Hit('contig_2_6', 'PF00001.1', 110, 300, 2, 190, 1, 195)
        rm = self._manager({'PF00001.1': [a, b, c]})
        rm.identifyAdjacentMarkerGenes()

        self.assertEqual(sorted(h.target_name for h in rm.markerHits['PF00001.1']),
                         ['contig_1_5', 'contig_1_7', 'contig_2_6'])

    def testMultipleMerges(self):
        """Several pairs of adjacent hits are merged; other markers are untouched."""
        a = Hit('contig_1_5', 'PF00001.1', 1, 120, 3, 118, 1, 120)
        b = Hit('contig_1_6', 'PF00001.1', 110, 300, 2, 190, 1, 195)
        c = Hit('contig_2_1', 'PF00001.1', 1, 100, 1, 100, 1, 100)
        d = Hit('contig_2_2', 'PF00001.1', 90, 250, 5, 160, 5, 165)
        e = Hit('contig_1_6', 'TIGR00001', 1, 100, 1, 100, 1, 100)
        rm = self._manager({'PF00001.1': [a, b, c, d], 'TIGR00001': [e]})
        rm.identifyAdjacentMarkerGenes()

        merged = {h.target_name: h for h in rm.markerHits['PF00001.1']}
        self.assertEqual(sorted(merged), ['contig_1_5&&contig_1_6', 'contig_2_1&&contig_2_2'])
        h = merged['contig_2_1&&contig_2_2']
        self.assertEqual((h.hmm_from, h.hmm_to, h.ali_from, h.ali_to, h.env_from, h.env_to),
                         (1, 250, 1, 160, 1, 165))
        self.assertEqual(rm.markerHits['TIGR00001'], [e])

    def testNonProdigalNames(self):
        """Gene ids without a numeric suffix are left as they are."""
        a = Hit('geneA_x', 'PF00001.1', 1, 120, 3, 118, 1, 120)
        b = Hit('geneA_y', 'PF00001.1', 110, 300, 2, 190, 1, 195)
        rm = self._manager({'PF00001.1': [a, b], 'PF00002.1': []})
        rm.identifyAdjacentMarkerGenes()
        self.assertEqual(rm.markerHits['PF00001.1'], [a, b])
        self.assertEqual(rm.markerHits['PF00002.1'], [])


if __name__ == "__main__":
    unittest.main()
