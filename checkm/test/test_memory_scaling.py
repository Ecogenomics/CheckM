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

"""Tests of the changes that keep memory use independent of the number of bins.

Bins are processed one at a time (HMM hits are parsed, used and discarded), and
bins share the marker sets and HMM models they have in common. These tests verify
that this gives the same results as holding the results of all bins in memory.
"""

import gzip
import os
import shutil
import tempfile
import unittest

from checkm.defaultValues import DefaultValues
from checkm.hmmerAligner import HmmerAligner
from checkm.hmmerModelParser import HmmModel
from checkm.markerSets import MarkerSet, BinMarkerSets
from checkm.resultsParser import ResultsParser
from checkm.aminoAcidIdentity import AminoAcidIdentity

MARKERS = ['PF00001.1', 'PF00002.1', 'TIGR00003']


def haveReferenceData():
    return os.path.exists(DefaultValues.PFAM_CLAN_FILE)


def domtbloutLine(geneId, marker, score, eValue='1e-30'):
    """Line of a HMMER domain table (domtblout) for a full-length hit."""
    fields = [geneId, '-', '100', marker, marker, '100', eValue, str(score), '0.1', '1', '1',
              eValue, eValue, str(score), '0.1', '1', '100', '1', '100', '1', '100', '0.99',
              '# 1 # 300 # 1 # ID=1_1']
    return ' '.join(fields) + '\n'


class Hit():
    def __init__(self, target_name, full_e_value):
        self.target_name = target_name
        self.full_e_value = full_e_value


class VerifySharedMarkerSets(unittest.TestCase):
    LINE = "%s\t2\tUID1\tk__Bacteria\t100\t[{'PF00001.1', 'TIGR00398'}, {'PF00002.1'}, {'TIGR00399'}]" \
           "\tUID0\troot\t500\t[{'PF00001.1'}]\n"

    def _read(self, binId, sharedMarkerSets=None):
        bms = BinMarkerSets(binId, BinMarkerSets.TREE_MARKER_SET)
        bms.read(self.LINE % binId, sharedMarkerSets)
        return bms

    def testSharedMarkerSets(self):
        """Bins with the same marker sets share them only if read with the same dictionary."""
        shared = {}
        a = self._read('binA', shared)
        b = self._read('binB', shared)
        c = self._read('binC')

        self.assertEqual(len(shared), 2)
        for i in range(2):
            self.assertIs(a.markerSets[i], b.markerSets[i])
            self.assertIsNot(a.markerSets[i], c.markerSets[i])
            self.assertEqual(str(a.markerSets[i]), str(c.markerSets[i]))

        self.assertEqual(a.markerSets[0].UID, 'UID1')
        self.assertEqual(a.markerSets[0].numGenomes, 100)
        self.assertEqual(a.markerSets[1].lineageStr, 'root')

    def testRemoveMarkersFromSharedMarkerSets(self):
        """Removing markers for every bin gives the result of removing them once."""
        shared = {}
        binMarkerSets = [self._read(binId, shared) for binId in ['binA', 'binB', 'binC']]
        for bms in binMarkerSets:
            bms.removeMarkers(DefaultValues.MARKERS_TO_EXCLUDE)

        unshared = self._read('binD')
        unshared.removeMarkers(DefaultValues.MARKERS_TO_EXCLUDE)

        for bms in binMarkerSets:
            self.assertEqual(bms.markerSets[0].markerSet, [{'PF00001.1'}, {'PF00002.1'}])
            self.assertEqual(bms.markerSets[0].markerSet, unshared.markerSets[0].markerSet)
            self.assertEqual(bms.getMarkerGenes(), {'PF00001.1', 'PF00002.1'})

    def testRemoveMarkersKeepsUnchangedSets(self):
        """Sets without markers to remove are not copied."""
        unchanged = frozenset(['PF00001.1'])
        ms = MarkerSet('UID1', 'k__Bacteria', 100, [unchanged, frozenset(['PF00002.1', 'TIGR00398']), frozenset(['TIGR00399'])])
        ms.removeMarkers({'TIGR00398', 'TIGR00399'})

        self.assertEqual(ms.markerSet, [{'PF00001.1'}, {'PF00002.1'}])
        self.assertIs(ms.markerSet[0], unchanged)


class VerifyMultipleHitExtraction(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        self.bms = BinMarkerSets('bin1', BinMarkerSets.HMM_MODELS_SET)
        self.bms.addMarkerSet(MarkerSet(0, 'N/A', -1, [{'PF00001.1', 'PF00002.1'}]))

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testNoMultipleHits(self):
        """The genes of a bin are not read if no marker has multiple hits."""
        markerHits = {'PF00001.1': [Hit('c_1', 1e-10)],
                      'TIGR00003': [Hit('c_2', 1e-10), Hit('c_7', 1e-20)]}  # not in the marker set
        HA = HmmerAligner(1)

        # there is no gene file, so this fails if an attempt is made to read it
        self.assertEqual(len(HA._extractMarkersWithMultipleHits(self.tmpDir, 'bin1', markerHits, self.bms)), 0)

    def testMultipleHits(self):
        binDir = os.path.join(self.tmpDir, 'bins', 'bin1')
        os.makedirs(binDir)
        with gzip.open(os.path.join(binDir, DefaultValues.PRODIGAL_AA), 'wt') as f:
            f.write('>c_1 # 1 # 30\nMKA*\n>c_5 # 40 # 70\nMKC\n>c_9\nMKD*\n')

        markerHits = {'PF00001.1': [Hit('c_1', 1e-30), Hit('c_5', 1e-10)],
                      'PF00002.1': [Hit('c_9', 1e-10)]}
        HA = HmmerAligner(1)
        multipleHits = HA._extractMarkersWithMultipleHits(self.tmpDir, 'bin1', markerHits, self.bms)

        # sequences are ordered from the highest to the lowest e-value
        self.assertEqual(list(multipleHits), ['PF00001.1'])
        self.assertEqual(list(multipleHits['PF00001.1']['bin1'].items()), [('c_5', 'MKC'), ('c_1', 'MKA')])


@unittest.skipUnless(haveReferenceData(), 'CheckM reference data is not available')
class VerifyReportResults(unittest.TestCase):
    """Reporting bins one at a time gives the output of reporting bins held in memory."""

    NUM_BINS = 23

    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        os.makedirs(os.path.join(self.tmpDir, 'storage'))

        models = {m: HmmModel({'name': m, 'acc': m, 'leng': 100, 'ga': (25.0, 25.0)}) for m in MARKERS}

        self.binIds = ['bin_%02d' % i for i in range(self.NUM_BINS)]
        self.binIdToModels = {}
        self.binIdToBinMarkerSets = {}
        binStats = {}
        for i, binId in enumerate(self.binIds):
            binDir = os.path.join(self.tmpDir, 'bins', binId)
            os.makedirs(binDir)
            with gzip.open(os.path.join(binDir, DefaultValues.HMMER_TABLE_OUT), 'wt') as f:
                f.write('# target name\n')
                if i % 5 != 0:
                    f.write(domtbloutLine('c1_1', 'PF00001.1', 100))
                if i % 3 == 0:
                    # second copy of a marker
                    f.write(domtbloutLine('c2_%d' % (i + 5), 'PF00001.1', 80))
                if i % 2 == 0:
                    f.write(domtbloutLine('c1_9', 'PF00002.1', 60))
                if i % 4 == 0:
                    # adjacent genes are combined into a single hit
                    f.write(domtbloutLine('c3_4', 'TIGR00003', 50))
                    f.write(domtbloutLine('c3_5', 'TIGR00003', 55))
                f.write(domtbloutLine('c4_1', 'PF00002.1', 5))  # below the model threshold

            self.binIdToModels[binId] = models

            bms = BinMarkerSets(binId, BinMarkerSets.HMM_MODELS_SET)
            bms.addMarkerSet(MarkerSet(0, 'N/A', -1, [{'PF00001.1', 'PF00002.1'}, {'TIGR00003'}]))
            self.binIdToBinMarkerSets[binId] = bms

            binStats[binId] = {'Genome size': 1000 + i, '# ambiguous bases': 0, '# scaffolds': 2, '# contigs': 3,
                               'N50 (scaffolds)': 500, 'N50 (contigs)': 400, 'Mean scaffold length': 500.0,
                               'Mean contig length': 333.3, 'Longest scaffold': 600, 'Longest contig': 500,
                               'GC': 0.5, 'GC std': 0.01, 'Coding density': 0.9, 'Translation table': 11,
                               '# predicted genes': 10 + i}

        # bins are not in sorted order, as output must be ordered by bin ID regardless
        self.binIdToModels = {binId: self.binIdToModels[binId] for binId in reversed(self.binIds)}

        with gzip.open(os.path.join(self.tmpDir, 'storage', DefaultValues.BIN_STATS_OUT), 'wt') as f:
            for binId in self.binIds:
                f.write(binId + '\t' + str(binStats[binId]) + '\n')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def _cacheFiles(self):
        contents = []
        for name in [DefaultValues.BIN_STATS_EXT_OUT, DefaultValues.MARKER_GENE_STATS]:
            path = os.path.join(self.tmpDir, 'storage', name)
            with gzip.open(path, 'rt') as f:
                contents.append(f.read())
            os.remove(path)
        return contents

    def _inMemory(self, threads, outputFormat, bTabTable):
        outFile = os.path.join(self.tmpDir, 'in_memory.txt')
        RP = ResultsParser(self.binIdToModels, threads)
        RP.analyseResults(self.tmpDir, DefaultValues.BIN_STATS_OUT, DefaultValues.HMMER_TABLE_OUT)
        RP.printSummary(outputFormat, AminoAcidIdentity(), self.binIdToBinMarkerSets, False, None, bTabTable, outFile, self.tmpDir)
        RP.cacheResults(self.tmpDir, self.binIdToBinMarkerSets, False)
        with open(outFile) as f:
            return [f.read()] + self._cacheFiles()

    def _oneAtATime(self, threads, outputFormat, bTabTable):
        outFile = os.path.join(self.tmpDir, 'one_at_a_time.txt')
        RP = ResultsParser(self.binIdToModels, threads)
        RP.reportResults(self.tmpDir, DefaultValues.BIN_STATS_OUT, DefaultValues.HMMER_TABLE_OUT,
                         outputFormat, AminoAcidIdentity(), self.binIdToBinMarkerSets, False, None, bTabTable, outFile)
        self.assertEqual(RP.results, {})
        with open(outFile) as f:
            return [f.read()] + self._cacheFiles()

    def testSameOutput(self):
        expected = self._inMemory(1, 2, True)
        self.assertEqual(expected[0].count('\n'), self.NUM_BINS + 1)
        self.assertEqual(expected[1].count('\n'), self.NUM_BINS)
        self.assertIn("'GCN2': ['PF00001']", expected[1])

        # more bins than are parsed ahead with 2 or 4 threads
        for threads in [1, 2, 4]:
            self.assertEqual(self._oneAtATime(threads, 2, True), expected)
            self.assertEqual(self._inMemory(threads, 2, True), expected)

    def testOtherOutputFormats(self):
        for outputFormat, bTabTable in [(1, False), (2, False), (3, True), (5, True), (6, True), (8, True)]:
            expected = self._inMemory(1, outputFormat, bTabTable)
            self.assertTrue(expected[0])
            self.assertEqual(self._oneAtATime(3, outputFormat, bTabTable), expected)

    def testBinsYieldedInOrder(self):
        RP = ResultsParser(self.binIdToModels, 4)

        class Progress():
            count = 0

            def update(self, n):
                self.count += n

        progress = Progress()
        args = (False, DefaultValues.E_VAL, DefaultValues.LENGTH, False, None)
        binIds = sorted(self.binIds)
        results = list(RP._iterBinResults(self.tmpDir, DefaultValues.HMMER_TABLE_OUT, binIds, False, args, progress))

        self.assertEqual([rm.binId for rm in results], binIds)
        self.assertEqual(progress.count, self.NUM_BINS)

        # bin_00 has only the second copy of PF00001.1 and adjacent hits to TIGR00003
        self.assertEqual([h.target_name for h in results[0].markerHits['PF00001.1']], ['c2_5'])
        self.assertEqual([h.target_name for h in results[0].markerHits['TIGR00003']], ['c3_4&&c3_5'])
        self.assertEqual(sorted(h.target_name for h in results[3].markerHits['PF00001.1']), ['c1_1', 'c2_8'])
        self.assertNotIn('PF00002.1', results[3].markerHits)


class VerifyDeepTree(unittest.TestCase):
    def testDeeplyNestedTree(self):
        """Trees more deeply nested than Python's default recursion limit can be read."""
        import sys
        from checkm.treeParser import readTree

        depth = 3000
        newick = '(' * depth + 'bin_0:0.1' + ''.join(',bin_%d:0.1):0.1' % i for i in range(1, depth + 1)) + ';\n'

        tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        origLimit = sys.getrecursionlimit()
        try:
            treeFile = os.path.join(tmpDir, 'deep.tre.gz')
            with gzip.open(treeFile, 'wt') as f:
                f.write(newick)

            tree = readTree(treeFile)
            self.assertEqual(len(tree.leaf_nodes()), depth + 1)
            self.assertEqual(tree.find_node_with_taxon_label('bin_0').level(), depth)
        finally:
            sys.setrecursionlimit(origLimit)
            shutil.rmtree(tmpDir)


if __name__ == '__main__':
    unittest.main()
