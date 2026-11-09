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
import logging
import unittest

from checkm.markerSets import MarkerSet, BinMarkerSets


class VerifyGenomeCheck(unittest.TestCase):
    def setUp(self):
        logging.getLogger('timestamp').disabled = True

    def tearDown(self):
        logging.getLogger('timestamp').disabled = False

    def testCollocatedSets(self):
        ms = MarkerSet('UID1', 'k__Bacteria', 10, [set(['a', 'b']), set(['c'])])
        hits = {'a': ['h1', 'h2'], 'c': ['h3']}
        comp, cont = ms.genomeCheck(hits, False)
        self.assertAlmostEqual(comp, 100 * (0.5 + 1.0) / 2)
        self.assertAlmostEqual(cont, 100 * (0.5 + 0.0) / 2)

    def testIndividualMarkers(self):
        ms = MarkerSet('UID1', 'k__Bacteria', 10, [set(['a', 'b']), set(['c', 'd'])])
        hits = {'a': ['h1', 'h2', 'h3'], 'c': ['h4']}
        comp, cont = ms.genomeCheck(hits, True)
        self.assertAlmostEqual(comp, 50.0)
        self.assertAlmostEqual(cont, 50.0)

    def testIndividualMarkersEmptyHitLists(self):
        """Empty hit lists are absent markers, not markers with -1 extra copies."""
        ms = MarkerSet('UID1', 'k__Bacteria', 10, [set(['a', 'b']), set(['c', 'd'])])
        hits = {'a': ['h1', 'h2', 'h3'], 'b': [], 'c': ['h4'], 'd': []}
        self.assertEqual(ms.genomeCheck(hits, True), (50.0, 50.0))
        self.assertEqual(ms.genomeCheck(hits, False), (50.0, 50.0))

    def testEmptyMarkerSet(self):
        """All markers removed (e.g., qa --exclude_markers) must not divide by zero."""
        ms = MarkerSet('UID1', 'k__Bacteria', 10, [set(['a', 'b']), set(['c'])])
        ms.removeMarkers(set(['a', 'b', 'c']))
        self.assertEqual(ms.markerSet, [])
        self.assertEqual(ms.genomeCheck({'a': ['h1']}, False), (0.0, 0.0))
        self.assertEqual(ms.genomeCheck({'a': ['h1']}, True), (0.0, 0.0))

        ms = MarkerSet(0, 'N/A', -1, [set()])
        self.assertEqual(ms.genomeCheck({}, False), (0.0, 0.0))
        self.assertEqual(ms.genomeCheck({}, True), (0.0, 0.0))

    def testEmptyCollocatedSetIgnored(self):
        ms = MarkerSet('UID1', 'k__Bacteria', 10, [set(['a', 'b']), set()])
        self.assertEqual(ms.genomeCheck({'a': ['h1']}, False), (50.0, 0.0))


class VerifySelectedMarkerSet(unittest.TestCase):
    def setUp(self):
        logging.getLogger('timestamp').disabled = True

    def tearDown(self):
        logging.getLogger('timestamp').disabled = False

    def _binMarkerSets(self, uids):
        bms = BinMarkerSets('bin1', BinMarkerSets.TREE_MARKER_SET)
        for uid in uids:
            bms.addMarkerSet(MarkerSet(uid, uid, 10, [set([uid + '_marker'])]))
        return bms

    def testFullTree(self):
        """The marker set selected for the most specific node is used."""
        bms = self._binMarkerSets(['UID5', 'UID4', 'UID3', 'UID1'])
        bms.setLineageSpecificSelectedMarkerSet({'UID5': 'UID4', 'UID4': 'UID4', 'UID3': 'UID3', 'UID1': 'UID1'})
        self.assertEqual(bms.selectedMarkerSet().UID, 'UID4')

    def testReducedTreeFollowChain(self):
        """A selected set missing from the reduced tree is replaced by the set it selects."""
        bms = self._binMarkerSets(['UID5', 'UID3', 'UID1'])
        bms.setLineageSpecificSelectedMarkerSet({'UID5': 'UID4', 'UID4': 'UID3', 'UID3': 'UID3', 'UID1': 'UID1'})
        self.assertEqual(bms.selectedMarkerSet().UID, 'UID3')

    def testReducedTreeUID3097(self):
        """UID3097 selects UID3076, which selects itself and is not in the reduced tree.

        This previously looped forever. The next (less specific) marker set of
        the bin is now tried instead.
        """
        bms = self._binMarkerSets(['UID3097', 'UID2', 'UID1'])
        selectedMarkerSetMap = {'UID3097': 'UID3076', 'UID3076': 'UID3076', 'UID2': 'UID2', 'UID1': 'UID1'}
        bms.setLineageSpecificSelectedMarkerSet(selectedMarkerSetMap)
        self.assertEqual(bms.selectedMarkerSet().UID, 'UID2')

    def testReducedTreeLongerCycle(self):
        bms = self._binMarkerSets(['UID9', 'UID7', 'UID1'])
        selectedMarkerSetMap = {'UID9': 'UID8', 'UID8': 'UID6', 'UID6': 'UID8',
                                'UID7': 'UID6', 'UID1': 'UID1'}
        bms.setLineageSpecificSelectedMarkerSet(selectedMarkerSetMap)
        self.assertEqual(bms.selectedMarkerSet().UID, 'UID1')

    def testNoSelectableMarkerSet(self):
        bms = self._binMarkerSets(['UID3097', 'UID2'])
        selectedMarkerSetMap = {'UID3097': 'UID3076', 'UID3076': 'UID3076', 'UID2': 'UID3076'}
        with self.assertRaises(SystemExit) as cm:
            bms.setLineageSpecificSelectedMarkerSet(selectedMarkerSetMap)
        self.assertEqual(cm.exception.code, 1)

    def testSelectedMarkerSetData(self):
        """Every node in the reference data selects a marker set in a finite number of steps."""
        from checkm.defaultValues import DefaultValues
        from checkm.markerSets import MarkerSetParser

        if not os.path.exists(DefaultValues.SELECTED_MARKER_SETS):
            self.skipTest('CheckM reference data not available')

        selectedMarkerSetMap = MarkerSetParser().parseSelectedMarkerSetMap()
        self.assertEqual(selectedMarkerSetMap['UID3097'], 'UID3076')
        self.assertEqual(selectedMarkerSetMap['UID3076'], 'UID3076')

        bms = self._binMarkerSets(['UID3097', 'UID2', 'UID1'])
        bms.setLineageSpecificSelectedMarkerSet(selectedMarkerSetMap)
        self.assertIn(bms.selectedMarkerSet().UID, ['UID2', 'UID1'])


class VerifyMarkerSetFile(unittest.TestCase):
    def testWriteRead(self):
        """Marker sets are written with sorted markers and can be read back."""
        markers = [set(['TIGR00002', 'PF00001.1', 'b']), set(['a']), set()]
        ms = MarkerSet('UID2', 'k__Bacteria', 10, markers)
        self.assertEqual(str(ms), "UID2\tk__Bacteria\t10\t[{'PF00001.1', 'TIGR00002', 'b'}, {'a'}, set()]")

        bms = BinMarkerSets('bin1', BinMarkerSets.TREE_MARKER_SET)
        bms.addMarkerSet(ms)
        bms.addMarkerSet(MarkerSet('UID1', 'root', 20, [set(['x', 'y'])]))

        out = io.StringIO()
        bms.write(out)
        line = out.getvalue()
        self.assertEqual(line, "bin1\t2\tUID2\tk__Bacteria\t10\t[{'PF00001.1', 'TIGR00002', 'b'}, {'a'}, set()]"
                               "\tUID1\troot\t20\t[{'x', 'y'}]\n")

        readBms = BinMarkerSets('bin1', BinMarkerSets.TREE_MARKER_SET)
        readBms.read(line)
        self.assertEqual([m.UID for m in readBms.markerSets], ['UID2', 'UID1'])
        self.assertEqual(readBms.markerSets[0].markerSet, markers)
        self.assertEqual(readBms.markerSets[1].markerSet, [set(['x', 'y'])])
        self.assertEqual(readBms.markerSets[1].numGenomes, 20)


if __name__ == "__main__":
    unittest.main()
