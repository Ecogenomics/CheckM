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

import ast
import unittest

from checkm.markerSets import MarkerSet, parseMarkerSetStr, parseSetStr


class VerifyMarkerSetParsing(unittest.TestCase):
    def testFormats(self):
        """Verify marker set strings in all formats written by CheckM are parsed like Python literals."""

        strs = ["[{'PF00001.1', 'TIGR00001'}, {'PF00002.3'}]",            # CheckM >= 1.1 (Python 3)
                "[set(['PF00001.1', 'TIGR00001']), set(['PF00002.3'])]",  # CheckM < 1.1 (Python 2) and reference data
                '[]', '[set()]', '[set([])]', "[{'a'}, set(), set(['b', 'c'])]",
                "[{'a' , 'b' }]", '[{"a"}]', '[{"x\'y"}]', "[{'a,b', 'c'}]", "[{'a\\\\b'}]"]

        for s in strs:
            expected = [set(x) for x in eval(s.replace('set([])', 'set()'))]
            self.assertEqual(parseMarkerSetStr(s), expected, s)

    def testSingleSets(self):
        """Verify single set strings are parsed."""
        self.assertEqual(parseSetStr('set([])'), set())
        self.assertEqual(parseSetStr("set(['pfam08459'])"), {'pfam08459'})
        self.assertEqual(parseSetStr("set(['a', 'b'])"), {'a', 'b'})
        self.assertEqual(parseSetStr("{'a', 'b'}"), {'a', 'b'})

    def testRoundTrip(self):
        """Verify marker sets written by MarkerSet are parsed back unchanged."""
        markers = [{'PF00001.1', 'TIGR00001'}, set(), {'PF00002.3'}]
        ms = MarkerSet('UID1', 'k__Bacteria', 10, markers)
        markerSetStr = str(ms).split('\t')[3]
        self.assertEqual(parseMarkerSetStr(markerSetStr), markers)

    def testNoCodeEvaluation(self):
        """Verify file contents are not evaluated as code."""
        with self.assertRaises((ValueError, SyntaxError)):
            parseMarkerSetStr("[__import__('os').getcwd()]")

    def testSharedSetsAreImmutable(self):
        """Verify parsed sets can't be modified, as they may be shared between marker sets."""
        sets = parseMarkerSetStr("[{'a', 'b'}]")
        self.assertIsInstance(sets[0], frozenset)
        self.assertIsNot(parseMarkerSetStr("[{'a', 'b'}]"), sets)


if __name__ == "__main__":
    unittest.main()
