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
from types import SimpleNamespace

from checkm.main import OptionsParser


class VerifyPplacerThreads(unittest.TestCase):
    def testDefault(self):
        """pplacer threads default to the number of threads, to a maximum of 64."""
        op = OptionsParser()
        self.assertEqual(op._pplacerThreads(SimpleNamespace(threads=1, pplacer_threads=None)), 1)
        self.assertEqual(op._pplacerThreads(SimpleNamespace(threads=48, pplacer_threads=None)), 48)
        self.assertEqual(op._pplacerThreads(SimpleNamespace(threads=96, pplacer_threads=None)), 64)
        self.assertEqual(op._pplacerThreads(SimpleNamespace(threads=8)), 8)

    def testExplicit(self):
        """An explicit number of pplacer threads is used, to a maximum of 64."""
        op = OptionsParser()
        self.assertEqual(op._pplacerThreads(SimpleNamespace(threads=64, pplacer_threads=4)), 4)
        with self.assertLogs('timestamp', level='WARNING'):
            self.assertEqual(op._pplacerThreads(SimpleNamespace(threads=8, pplacer_threads=128)), 64)


if __name__ == "__main__":
    unittest.main()
