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
