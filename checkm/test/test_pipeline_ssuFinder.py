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

from checkm.ssuFinder import SSU_Finder


NHMMER_OUTPUT = """# nhmmer :: search a DNA model, alignment, or sequence against a DNA database
>> contig_1
    score  bias    Evalue   hmmfrom    hmm to     alifrom    ali to      envfrom    env to       sq len      acc
   ------ ----- ---------   -------   -------    --------- ---------    --------- ---------    ---------    ----
 ! 1335.0  27.2         0        29      1514 ..      2000       501 ..      2010       490 ..      5000    0.93

>> contig_10
    score  bias    Evalue   hmmfrom    hmm to     alifrom    ali to      envfrom    env to       sq len      acc
   ------ ----- ---------   -------   -------    --------- ---------    --------- ---------    ---------    ----
 !  520.0  42.1  8.9e-155        35       775 ..       101       900 ..        90       910 ..      6006    0.92
 !   20.0   2.1     1e-03         1        50 ..      1001      1050 ..       990      1060 ..      6006    0.92

"""


def hit(start, end, bRev, domain='bacteria', evalue='0'):
    return [domain, evalue, str(start), str(end), str(end - start + 1), str(bRev)]


class VerifySSUFinder(unittest.TestCase):
    def setUp(self):
        self.ssuFinder = SSU_Finder(threads=1)

    def testReadHits(self):
        """Verify parsing of nhmmer hits, including strand and inclusive lengths."""
        tmpDir = tempfile.mkdtemp()
        try:
            resultsFile = os.path.join(tmpDir, 'ssu.bacteria.txt')
            with open(resultsFile, 'w') as f:
                f.write(NHMMER_OUTPUT)

            seqInfo = self.ssuFinder._readHits(resultsFile, 'bacteria', 1e-5)
        finally:
            shutil.rmtree(tmpDir)

        # reverse strand hit: coordinates are swapped so start < end
        self.assertEqual(seqInfo['contig_1'], [['bacteria', '0', '501', '2000', '1500', 'True']])

        # forward strand hit; second hit is above the e-value threshold
        self.assertEqual(seqInfo['contig_10'], [['bacteria', '8.9e-155', '101', '900', '800', 'False']])

    def testAddHitConcatenatesSameStrand(self):
        """Verify nearby hits on the same strand are concatenated."""
        hits = {}
        self.ssuFinder._addHit(hits, 'contig_1', hit(100, 500, False), 200)
        self.ssuFinder._addHit(hits, 'contig_1', hit(600, 1500, False), 200)

        self.assertEqual(hits, {'contig_1': hit(100, 1500, False)})
        self.assertEqual(hits['contig_1'][4], '1401')

    def testAddHitDoesNotConcatenateOppositeStrands(self):
        """Verify nearby hits on opposite strands are kept as separate hits."""
        hits = {}
        self.ssuFinder._addHit(hits, 'contig_1', hit(100, 500, False), 200)
        self.ssuFinder._addHit(hits, 'contig_1', hit(600, 1500, True), 200)

        self.assertEqual(hits, {'contig_1': hit(100, 500, False),
                                'contig_1-#2': hit(600, 1500, True)})

        # same result with the strands swapped
        hits = {}
        self.ssuFinder._addHit(hits, 'contig_1', hit(100, 500, True), 200)
        self.ssuFinder._addHit(hits, 'contig_1', hit(600, 1500, False), 200)

        self.assertEqual(hits, {'contig_1': hit(100, 500, True),
                                'contig_1-#2': hit(600, 1500, False)})

    def testAddDomainHitExactContigMatch(self):
        """Verify hits are only compared with hits on the same contig (contig_1 is not contig_10)."""
        bestHits = {}
        self.ssuFinder._addDomainHit(bestHits, 'contig_10', hit(100, 1600, False, 'bacteria'))

        # overlapping, shorter hit on a different contig whose ID is a prefix of contig_10
        self.ssuFinder._addDomainHit(bestHits, 'contig_1', hit(200, 1000, False, 'archaea'))

        self.assertEqual(bestHits, {'contig_10-#1': hit(100, 1600, False, 'bacteria'),
                                    'contig_1-#1': hit(200, 1000, False, 'archaea')})

        # overlapping, longer hit on contig_1 replaces the hit on contig_1 only
        self.ssuFinder._addDomainHit(bestHits, 'contig_1', hit(150, 1200, False, 'euk'))
        self.assertEqual(bestHits, {'contig_10-#1': hit(100, 1600, False, 'bacteria'),
                                    'contig_1-#1': hit(150, 1200, False, 'euk')})

        # non-overlapping hit on contig_1 is added as a new hit
        self.ssuFinder._addDomainHit(bestHits, 'contig_1', hit(3000, 4000, True, 'bacteria'))
        self.assertEqual(bestHits['contig_1-#2'], hit(3000, 4000, True, 'bacteria'))
        self.assertEqual(len(bestHits), 3)

    def testExtractSeq(self):
        """Verify extraction uses 1-based inclusive coordinates and reverse complements reverse strand hits."""
        seq = 'AACCGGTTAC'

        self.assertEqual(self.ssuFinder._ssuSeq(seq, 1, 10, False), seq)
        self.assertEqual(self.ssuFinder._ssuSeq(seq, 3, 6, False), 'CCGG')
        self.assertEqual(self.ssuFinder._ssuSeq(seq, 7, 10, False), 'TTAC')
        self.assertEqual(self.ssuFinder._ssuSeq(seq, 7, 10, True), 'GTAA')
        self.assertEqual(self.ssuFinder._ssuSeq('ACGTNacgtn', 1, 10, True), 'nacgtNACGT')


if __name__ == "__main__":
    unittest.main()
