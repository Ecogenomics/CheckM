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

"""Tests for the order in which bins are processed by MarkerGeneFinder (memory budget for gene calling)."""

import gzip
import heapq
import os
import random
import shutil
import tempfile
import unittest

from checkm.defaultValues import DefaultValues
from checkm.markerGeneFinder import (_ReadyBins, _genomeSize, _geneCallingMemoryBudget,
                                     _needsGeneCalling, _prodigalMemory)


class VerifyGenomeSize(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def testPlainAndGzip(self):
        data = b'>seq\n' + b'ACGT' * 25000 + b'\n'
        plain = os.path.join(self.tmpDir, 'g.fna')
        with open(plain, 'wb') as f:
            f.write(data)
        compressed = os.path.join(self.tmpDir, 'g.fna.gz')
        with gzip.open(compressed, 'wb') as f:
            f.write(data)

        self.assertEqual(_genomeSize(plain), len(data))
        self.assertEqual(_genomeSize(compressed), len(data))   # uncompressed size, not the file size
        self.assertEqual(_genomeSize(os.path.join(self.tmpDir, 'missing.fna')), 0)

    def testNeedsGeneCalling(self):
        outDir = self.tmpDir
        binFile = os.path.join(self.tmpDir, 'bin_1.fna')
        binDir = os.path.join(outDir, 'bins', 'bin_1')
        os.makedirs(binDir)

        self.assertTrue(_needsGeneCalling(binFile, outDir, False, False))
        self.assertFalse(_needsGeneCalling(binFile, outDir, False, True))   # genes provided (--genes)

        with gzip.open(os.path.join(binDir, DefaultValues.PRODIGAL_AA), 'wt') as f:
            f.write('>gene_1\nMKV\n')
        self.assertFalse(_needsGeneCalling(binFile, outDir, False, False))
        self.assertTrue(_needsGeneCalling(binFile, outDir, True, False))    # nucleotide genes also required


class VerifyMemoryBudget(unittest.TestCase):
    def testBudget(self):
        self.assertEqual(_geneCallingMemoryBudget([], 64), 0)

        sizes = [1000000 * (i % 10 + 1) for i in range(500)]
        budget = _geneCallingMemoryBudget(sizes, 64)
        self.assertGreaterEqual(budget, _prodigalMemory(max(sizes)))
        self.assertLess(budget, sum(sorted(map(_prodigalMemory, sizes))[-64:]))   # below the 64 largest at once

        # a single thread: the budget allows the largest genome
        self.assertGreaterEqual(_geneCallingMemoryBudget(sizes, 1), _prodigalMemory(max(sizes)))

        # few genomes per thread: the run can't be shorter than the largest genome
        self.assertGreaterEqual(_geneCallingMemoryBudget([10, 20, 30], 64), _prodigalMemory(30))


class VerifyReadyBins(unittest.TestCase):
    def testNoGeneCallingKeepsOrder(self):
        ready = _ReadyBins(0)
        ready.add(3, ['c1', 'c2'], {'c1': 0, 'c2': 0})
        ready.add(1, ['a1', 'a2', 'a3'], {'a1': 0, 'a2': 0, 'a3': 0})

        order = []
        while ready:
            groupIdx, binFile, memory = ready.pop()
            self.assertEqual(memory, 0)
            order.append((groupIdx, binFile))

        # groups in the order they became ready, bins in the given (largest first) order
        self.assertEqual(order, [(3, 'c1'), (3, 'c2'), (1, 'a1'), (1, 'a2'), (1, 'a3')])
        self.assertIsNone(ready.pop())

    def testLargestBinThatFits(self):
        memory = {'b100': 100, 'b60': 60, 'b50': 50, 'b10': 10, 'g0': 0}
        ready = _ReadyBins(120)
        ready.add(0, ['b100', 'b60', 'b50', 'b10', 'g0'], memory)

        self.assertEqual(ready.pop()[1], 'b100')
        self.assertEqual(ready.pop()[1], 'b10')     # 20 left: the largest bin that fits
        self.assertEqual(ready.pop()[1], 'g0')      # genes already called: no memory needed
        self.assertIsNone(ready.pop())              # nothing fits: wait
        self.assertEqual(ready.memoryInUse, 110)

        ready.release(100)
        self.assertEqual(ready.pop()[1], 'b60')
        self.assertEqual(ready.pop()[1], 'b50')
        self.assertFalse(ready)

    def testEqualMemoryKeepsOrder(self):
        ready = _ReadyBins(1000)
        ready.add(0, ['x', 'y', 'z'], {'x': 10, 'y': 10, 'z': 10})
        self.assertEqual([ready.pop()[1] for _ in range(3)], ['x', 'y', 'z'])

    def testSimulatedRunStaysWithinBudget(self):
        """Process random genomes with a fixed number of threads, as MarkerGeneFinder does."""
        rng = random.Random(1)
        threads = 64
        sizes = {f'bin_{i}': int(rng.lognormvariate(15, 0.6)) for i in range(2000)}
        binFiles = sorted(sizes, key=lambda f: -sizes[f])
        memory = {f: _prodigalMemory(s) for f, s in sizes.items()}

        budget = _geneCallingMemoryBudget(list(sizes.values()), threads)
        ready = _ReadyBins(budget)
        ready.add(0, binFiles, memory)

        t = 0
        running = []
        started = []
        peak = 0
        while ready or running:
            while len(running) < threads:
                nextBin = ready.pop() if ready else None
                if not nextBin:
                    break
                _, binFile, mem = nextBin
                heapq.heappush(running, (t + sizes[binFile], binFile, mem))
                started.append(binFile)
            peak = max(peak, ready.memoryInUse)
            self.assertTrue(running)   # never stalls
            t, binFile, mem = heapq.heappop(running)
            ready.release(mem)

        self.assertEqual(sorted(started), sorted(binFiles))
        self.assertLessEqual(peak, budget)
        self.assertEqual(ready.memoryInUse, 0)

        # largest first without a budget would run the 64 largest genomes at once
        self.assertLess(peak, sum(memory[f] for f in binFiles[0:threads]))

        # the run takes little longer than the shortest possible run
        makespan = max(sum(sizes.values()) / threads, max(sizes.values()))
        self.assertLess(t, 1.15 * makespan)


class VerifyFrozenHeap(unittest.TestCase):
    def testWorkersStartedWithFrozenHeap(self):
        """Worker processes are started with the heap frozen, so they do not copy the memory of the parent."""
        import gc
        from checkm.markerGeneFinder import MarkerGeneFinder
        freezeCounts = []

        def processGroups(*args, **kwargs):
            freezeCounts.append(gc.get_freeze_count())
            return {}, []

        mgf = MarkerGeneFinder(threads=2)
        mgf._groupBinsByHmms = lambda binFiles, markerFile: [(None, binFiles)]
        mgf._processGroups = processGroups

        tmpDir = tempfile.mkdtemp()
        try:
            gc.collect()
            initialCount = gc.get_freeze_count()   # Python may freeze some objects at startup
            mgf.find([os.path.join(tmpDir, 'bin_1.faa')], tmpDir, 'table.txt', 'hmmer.txt',
                     'markers.hmm', False, False, True)
        finally:
            shutil.rmtree(tmpDir)

        self.assertEqual(len(freezeCounts), 1)
        self.assertGreater(freezeCounts[0], initialCount)        # frozen while bins are processed
        self.assertLessEqual(gc.get_freeze_count(), initialCount)  # and unfrozen afterwards


if __name__ == '__main__':
    unittest.main()
