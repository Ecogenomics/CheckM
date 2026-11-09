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

import io
import os
import shutil
import logging
import tempfile
import unittest
from contextlib import redirect_stdout

from checkm.resultsParser import ResultsParser, ResultsManager
from checkm.markerSets import MarkerSet, BinMarkerSets
from checkm.aminoAcidIdentity import AminoAcidIdentity


class Hit(object):
    """Minimal stand-in for checkm.hmmer.HmmerHitDOM."""

    def __init__(self, target_name, query_accession,
                 hmm_from, hmm_to, ali_from, ali_to, env_from, env_to,
                 target_length=200, full_e_value=1e-30, full_score=100.0):
        self.target_name = target_name
        self.query_accession = query_accession
        self.target_length = target_length
        self.full_e_value = full_e_value
        self.full_score = full_score
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

    def testFragmentRuns(self):
        """A run of consecutive ORFs is combined into one hit, regardless of hit order."""
        import itertools

        hits = [Hit('contig_1_%d' % n, 'PF00001.1', 10 * n, 10 * n + 50, 1, 50, 1, 55,
                    target_length=60, full_e_value=1e-10 * n) for n in [28, 29, 30, 31]]
        for perm in itertools.permutations(hits):
            rm = self._manager({'PF00001.1': list(perm)})
            rm.identifyAdjacentMarkerGenes()

            combined = rm.markerHits['PF00001.1']
            self.assertEqual(len(combined), 1)
            h = combined[0]
            self.assertEqual(h.target_name, 'contig_1_28&&contig_1_29&&contig_1_30&&contig_1_31')
            self.assertEqual(h.target_length, 240)
            self.assertEqual((h.hmm_from, h.hmm_to), (280, 360))
            self.assertEqual([c.target_name for c in h.component_hits], ['contig_1_%d' % n for n in [28, 29, 30, 31]])

    def testFragmentRunGaps(self):
        """Runs are broken by gaps in ORF numbering; ORF numbers are compared numerically."""
        hits = [Hit('contig_1_%d' % n, 'PF00001.1', 1, 50, 1, 50, 1, 50) for n in [9, 10, 12, 13, 14, 20]]
        rm = self._manager({'PF00001.1': list(reversed(hits))})
        rm.identifyAdjacentMarkerGenes()

        self.assertEqual([h.target_name for h in rm.markerHits['PF00001.1']],
                         ['contig_1_9&&contig_1_10', 'contig_1_12&&contig_1_13&&contig_1_14', 'contig_1_20'])

    def testCombinedHitUsesBestHit(self):
        """Scores and e-values of a combined hit are those of the best hit."""
        a = Hit('contig_1_5', 'PF00001.1', 1, 120, 3, 118, 1, 120, full_e_value=1e-5, full_score=20.0)
        b = Hit('contig_1_6', 'PF00001.1', 110, 300, 2, 190, 1, 195, full_e_value=1e-50, full_score=200.0)
        rm = self._manager({'PF00001.1': [a, b]})
        rm.identifyAdjacentMarkerGenes()

        h = rm.markerHits['PF00001.1'][0]
        self.assertEqual((h.full_e_value, h.full_score), (1e-50, 200.0))

    def testNonProdigalNames(self):
        """Gene ids without a numeric suffix are left as they are."""
        a = Hit('geneA_x', 'PF00001.1', 1, 120, 3, 118, 1, 120)
        b = Hit('geneA_y', 'PF00001.1', 110, 300, 2, 190, 1, 195)
        rm = self._manager({'PF00001.1': [a, b], 'PF00002.1': []})
        rm.identifyAdjacentMarkerGenes()
        self.assertEqual(rm.markerHits['PF00001.1'], [a, b])
        self.assertEqual(rm.markerHits['PF00002.1'], [])


def hmmMarkerSets(binId, markers):
    """Bin marker sets as created for a marker file of HMMs."""
    bms = BinMarkerSets(binId, BinMarkerSets.HMM_MODELS_SET)
    bms.addMarkerSet(MarkerSet(0, 'N/A', -1, [set(markers)]))
    return bms


def capture(func, *args, **kwargs):
    out = io.StringIO()
    with redirect_stdout(out):
        rtn = func(*args, **kwargs)
    return out.getvalue(), rtn


class VerifyMarkerCounts(unittest.TestCase):
    def setUp(self):
        self.rm = ResultsManager('bin1', {})
        self.a = Hit('c_1', 'PF00001.1', 1, 100, 1, 100, 1, 100)
        self.b = Hit('c_5', 'PF00001.1', 1, 100, 1, 100, 1, 100)
        self.c = Hit('c_9', 'TIGR00002', 1, 100, 1, 100, 1, 100)
        self.rm.markerHits = {'PF00001.1': [self.a, self.b], 'TIGR00002': [self.c]}
        self.bms = hmmMarkerSets('bin1', ['PF00001.1', 'TIGR00002', 'TIGR00003', 'PF00004.1', 'A0001'])

    def testHitsToMarkerGeneDoesNotModifyHits(self):
        counts = self.rm.hitsToMarkerGene(self.bms.selectedMarkerSet())
        self.assertEqual(counts, {'PF00001.1': 2, 'TIGR00002': 1, 'TIGR00003': 0, 'PF00004.1': 0, 'A0001': 0})
        self.assertEqual(sorted(self.rm.markerHits), ['PF00001.1', 'TIGR00002'])

    def testIndividualMarkersAfterOutputFormat4(self):
        """qa -o 4 --individual_markers must not change completeness/contamination."""
        before = self.rm.geneCountsForSelectedMarkerSet(self.bms, True)
        capture(self.rm.printSummary, 4, None, self.bms, True)
        after = self.rm.geneCountsForSelectedMarkerSet(self.bms, True)
        self.assertEqual(before, after)
        self.assertEqual(after[0:6], [3, 1, 1, 0, 0, 0])
        self.assertAlmostEqual(after[6], 40.0)
        self.assertAlmostEqual(after[7], 20.0)

    def testOutputFormat4Sorted(self):
        out, _ = capture(self.rm.printSummary, 4, None, self.bms, False)
        lines = out.split('\n')
        self.assertEqual(lines[0], 'Node Id: 0; Marker lineage: N/A\tA0001\tPF00001.1\tPF00004.1\tTIGR00002\tTIGR00003')
        self.assertEqual(lines[1], 'bin1\t0\t2\t0\t1\t0')

    def testLineageOfHmmMarkerFile(self):
        """Marker sets from an HMM file have no UID, so none is appended to the lineage."""
        aai = AminoAcidIdentity()
        out, _ = capture(self.rm.printSummary, 1, aai, self.bms, False)
        self.assertTrue(out.startswith('bin1\tN/A\t-1\t5\t1\t'), out)

        # a lineage-specific set reports its UID
        bms = BinMarkerSets('bin1', BinMarkerSets.TAXONOMIC_MARKER_SET)
        bms.addMarkerSet(MarkerSet('UID1', 'k__Bacteria', 10, [set(['PF00001.1'])]))
        out, _ = capture(self.rm.printSummary, 1, aai, bms, False)
        self.assertTrue(out.startswith('bin1\tk__Bacteria (UID1)\t10\t1\t1\t'), out)


class VerifyOutputFormat9(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.tmpDir, 'bins', 'bin1'))
        self.genesFile = os.path.join(self.tmpDir, 'bins', 'bin1', 'genes.faa')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def _rows(self, rm, bms, table=None):
        out, _ = capture(rm.printSummary, 9, None, bms, False, None, table, self.tmpDir)
        return out

    def testMergedAndMultipleMarkers(self):
        with open(self.genesFile, 'w') as f:
            f.write('>c_1 # 1 # 300 # 1 # ID=1_1\nMAAA*\n')
            f.write('>c_2 # 400 # 900 # -1 # ID=1_2\nMCCCC*\n')
            f.write('>c_10 # 1000 # 1300 # 1 # ID=1_10\nMDD*\n')
            f.write('>c_3 # 950 # 990 # 1 # ID=1_3\nME*\n')

        # adjacent hits to c_1 and c_2 are combined; c_10 is hit by two markers
        a = Hit('c_1', 'PF00001.1', 1, 100, 5, 95, 1, 100)
        b = Hit('c_2', 'PF00001.1', 90, 250, 3, 160, 1, 165)
        c = Hit('c_10', 'PF00001.1', 1, 100, 2, 99, 1, 100)
        d = Hit('c_10', 'TIGR00002', 1, 100, 7, 80, 1, 100)
        e = Hit('c_3', 'TIGR00003', 1, 100, 7, 80, 1, 100)  # not in marker set

        rm = ResultsManager('bin1', {})
        rm.markerHits = {'TIGR00002': [d], 'PF00001.1': [a, b, c], 'TIGR00003': [e]}
        rm.identifyAdjacentMarkerGenes()
        bms = hmmMarkerSets('bin1', ['PF00001.1', 'TIGR00002'])

        rows = [r.split('\t') for r in self._rows(rm, bms).strip().split('\n')]
        self.assertEqual(rows, [
            ['bin1', 'c', '1', '1', '300', '1', '5', 'PF00001.1', '5', '95', 'MAAA*'],
            ['bin1', 'c', '2', '400', '900', '-1', '6', 'PF00001.1', '3', '160', 'MCCCC*'],
            ['bin1', 'c', '10', '1000', '1300', '1', '4', 'PF00001.1', '2', '99', 'MDD*'],
            ['bin1', 'c', '10', '1000', '1300', '1', '4', 'TIGR00002', '7', '80', 'MDD*'],
        ])

        fasta = self._rows(rm, bms, table=True).strip().split('\n')
        self.assertEqual(fasta[0], '>bin1 c geneId=1;start=1;end=300;strand=1;protlen=5 marker=PF00001.1;mstart=5;mend=95')
        self.assertEqual(fasta[1], 'MAAA*')
        self.assertEqual(len(fasta), 8)

    def testNonProdigalHeaders(self):
        """User-supplied protein files (--genes) need not use Prodigal headers."""
        with open(self.genesFile, 'w') as f:
            f.write('>geneA some description\nMAAA\n')
            f.write('>WP_000001.1 another protein\nMCC\n')
            f.write('>geneC\nMDD\n')

        rm = ResultsManager('bin1', {})
        rm.markerHits = {'PF00001.1': [Hit('geneA', 'PF00001.1', 1, 100, 5, 95, 1, 100)],
                         'TIGR00002': [Hit('WP_000001.1', 'TIGR00002', 1, 100, 2, 50, 1, 100)]}
        bms = hmmMarkerSets('bin1', ['PF00001.1', 'TIGR00002'])

        rows = [r.split('\t') for r in self._rows(rm, bms).strip().split('\n')]
        self.assertEqual(rows, [
            ['bin1', 'WP_000001.1', 'NA', 'NA', 'NA', 'NA', '3', 'TIGR00002', '2', '50', 'MCC'],
            ['bin1', 'geneA', 'NA', 'NA', 'NA', 'NA', '4', 'PF00001.1', '5', '95', 'MAAA'],
        ])


class VerifyParseHmmerResults(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        logging.getLogger('timestamp').disabled = True

    def tearDown(self):
        logging.getLogger('timestamp').disabled = False
        shutil.rmtree(self.tmpDir)

    def testMissingHmmerTable(self):
        """A missing HMMER table is an error, not a bin without marker genes."""
        rp = ResultsParser({'bin1': {}})
        rm = ResultsManager('bin1', {})
        missingFile = os.path.join(self.tmpDir, 'bins', 'bin1', 'hmmer.analyze.txt')
        with self.assertRaises(SystemExit) as cm:
            rp.parseHmmerResults(missingFile, rm, False)
        self.assertEqual(cm.exception.code, 1)

    def testMalformedHmmerTable(self):
        tableFile = os.path.join(self.tmpDir, 'hmmer.analyze.txt')
        with open(tableFile, 'w') as f:
            f.write('contig_1_1 - 100 PF00001 PF00001.1\n')

        rp = ResultsParser({'bin1': {}})
        rm = ResultsManager('bin1', {})
        with self.assertRaises(SystemExit) as cm:
            rp.parseHmmerResults(tableFile, rm, False)
        self.assertEqual(cm.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
