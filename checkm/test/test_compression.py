###############################################################################
#
# test_compression.py - test reading and writing gzip-compressed output files
#
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
import gzip
import shutil
import tarfile
import tempfile
import unittest

from checkm.util.compression import (resolvePath, openFile, compressFile, hasContent,
                                     isGzipped, compressedName, uncompressedName,
                                     archiveDirectory)
from checkm.util.seqUtils import readFasta, writeFasta
from checkm.aminoAcidIdentity import AminoAcidIdentity
from checkm.markerSets import MarkerSetParser, BinMarkerSets
from checkm.defaultValues import DefaultValues


class VerifyCompression(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm test ')

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def path(self, name):
        return os.path.join(self.tmpDir, name)

    def writePlain(self, name, text):
        with open(self.path(name), 'w') as f:
            f.write(text)
        return self.path(name)

    def testNames(self):
        self.assertEqual(compressedName('lineage.ms'), 'lineage.ms.gz')
        self.assertEqual(compressedName('lineage.ms.gz'), 'lineage.ms.gz')
        self.assertEqual(uncompressedName('genes.faa.gz'), 'genes.faa')
        self.assertEqual(uncompressedName('genes.faa'), 'genes.faa')

    def testResolvePath(self):
        # plain file given by its name
        plain = self.writePlain('a.txt', 'a\n')
        self.assertEqual(resolvePath(plain), plain)

        # compressed file given without '.gz' (e.g., lineage.ms -> lineage.ms.gz)
        with openFile(self.path('b.txt.gz'), 'wt') as f:
            f.write('b\n')
        self.assertEqual(resolvePath(self.path('b.txt')), self.path('b.txt.gz'))
        self.assertEqual(resolvePath(self.path('b.txt.gz')), self.path('b.txt.gz'))

        # uncompressed file of CheckM <= 1.2.x given by its new name (genes.faa.gz -> genes.faa)
        self.assertEqual(resolvePath(plain + '.gz'), plain)

        # the given path is preferred if both exist
        with openFile(plain + '.gz', 'wt') as f:
            f.write('a\n')
        self.assertEqual(resolvePath(plain), plain)
        self.assertEqual(resolvePath(plain + '.gz'), plain + '.gz')

        # missing files are returned unchanged
        self.assertEqual(resolvePath(self.path('missing.txt')), self.path('missing.txt'))
        self.assertEqual(resolvePath(self.path('missing.txt.gz')), self.path('missing.txt.gz'))

    def testOpenFileRoundTrip(self):
        text = 'line 1\nline 2\n'

        gzFile = self.path('x.txt.gz')
        with openFile(gzFile, 'wt') as f:
            f.write(text)
        self.assertTrue(isGzipped(gzFile))
        with gzip.open(gzFile, 'rt') as f:
            self.assertEqual(f.read(), text)
        with openFile(gzFile) as f:
            self.assertEqual(f.read(), text)
        with openFile(self.path('x.txt')) as f:  # resolved to x.txt.gz
            self.assertEqual(f.read(), text)

        # files without '.gz' are written uncompressed
        plainFile = self.path('y.txt')
        with openFile(plainFile, 'wt') as f:
            f.write(text)
        self.assertFalse(isGzipped(plainFile))
        with openFile(plainFile) as f:
            self.assertEqual(f.read(), text)

        # compression is determined from the content rather than the name
        misnamed = self.writePlain('z.txt.gz', text)
        with openFile(misnamed) as f:
            self.assertEqual(f.read(), text)

        # binary mode
        with openFile(gzFile, 'rb') as f:
            self.assertEqual(f.read(), text.encode())

    def testCompressFile(self):
        text = 'ACGT\n' * 1000
        src = self.writePlain('genes.faa.11', text)

        dst = compressFile(src, self.path('genes.faa.gz'), remove=False)
        self.assertEqual(dst, self.path('genes.faa.gz'))
        self.assertTrue(os.path.exists(src))
        self.assertTrue(isGzipped(dst))
        with gzip.open(dst, 'rt') as f:
            self.assertEqual(f.read(), text)

        # default destination, removing the source
        dst = compressFile(src)
        self.assertEqual(dst, src + '.gz')
        self.assertFalse(os.path.exists(src))
        with gzip.open(dst, 'rt') as f:
            self.assertEqual(f.read(), text)

        # a compressed source is not compressed twice
        dst2 = compressFile(dst, self.path('copy.gz'), remove=False)
        with gzip.open(dst2, 'rt') as f:
            self.assertEqual(f.read(), text)

    def testHasContent(self):
        emptyGz = self.path('empty.faa.gz')
        openFile(emptyGz, 'wt').close()
        self.assertGreater(os.path.getsize(emptyGz), 0)  # gzip header
        self.assertFalse(hasContent(emptyGz))
        self.assertFalse(hasContent(self.path('empty.faa')))  # resolved to empty.faa.gz

        with openFile(self.path('genes.faa.gz'), 'wt') as f:
            f.write('>g1\nMK\n')
        self.assertTrue(hasContent(self.path('genes.faa.gz')))

        self.assertFalse(hasContent(self.writePlain('empty.txt', '')))
        self.assertTrue(hasContent(self.writePlain('old.faa', '>g1\nMK\n')))
        self.assertTrue(hasContent(self.path('old.faa.gz')))  # resolved to old.faa
        self.assertFalse(hasContent(self.path('missing.faa.gz')))

    def testFasta(self):
        seqs = {'s1': 'ACGT', 's2': 'GGCC'}
        writeFasta(seqs, self.path('seqs.fna.gz'))
        self.assertTrue(isGzipped(self.path('seqs.fna.gz')))
        self.assertEqual(readFasta(self.path('seqs.fna.gz')), seqs)
        self.assertEqual(readFasta(self.path('seqs.fna')), seqs)

        writeFasta(seqs, self.path('plain.fna'))
        self.assertFalse(isGzipped(self.path('plain.fna')))
        self.assertEqual(readFasta(self.path('plain.fna.gz')), seqs)

    def testArchiveDirectory(self):
        inputDir = self.path('aai_qa')
        for binId, marker in [('bin_b', 'PF1'), ('bin_a', 'PF2'), ('bin_a', 'PF1')]:
            os.makedirs(os.path.join(inputDir, binId), exist_ok=True)
            with open(os.path.join(inputDir, binId, marker + '.masked.faa'), 'w') as f:
                f.write('>%s&&g1\nAC\n' % binId)

        archiveFile = self.path('aai_qa.tar.gz')
        self.writePlain('aai_qa.tar.gz', 'old archive')  # an existing archive is replaced
        archiveDirectory(inputDir, archiveFile)

        self.assertFalse(os.path.exists(inputDir))
        self.assertTrue(isGzipped(archiveFile))
        with tarfile.open(archiveFile, 'r:gz') as tar:
            self.assertEqual(tar.getnames(), ['bin_a/PF1.masked.faa', 'bin_a/PF2.masked.faa', 'bin_b/PF1.masked.faa'])
            self.assertEqual(tar.extractfile('bin_b/PF1.masked.faa').read(), b'>bin_b&&g1\nAC\n')

    def testMarkerFileType(self):
        markerFile = self.path('lineage.ms.gz')
        with openFile(markerFile, 'wt') as f:
            f.write(DefaultValues.LINEAGE_MARKER_FILE_HEADER + '\n')

        parser = MarkerSetParser()
        self.assertEqual(parser.markerFileType(markerFile), BinMarkerSets.TREE_MARKER_SET)
        self.assertEqual(parser.markerFileType(self.path('lineage.ms')), BinMarkerSets.TREE_MARKER_SET)


class VerifyAminoAcidIdentityArchive(unittest.TestCase):
    """Verify AAI is calculated from storage/aai_qa.tar.gz and from the storage/aai_qa directory of CheckM <= 1.2.x."""

    def setUp(self):
        self.outDir = tempfile.mkdtemp(prefix='checkm test ')

        alignments = {'bin_1': {'PF00001.masked.faa': {'bin_1&&g1': 'ACDEFGHIKL',
                                                       'bin_1&&g2': 'ACDEFGHIKM',
                                                       'bin_1&&g3': '--DEFGWIKL'},
                                'PF00002.masked.faa': {'bin_1&&g4': 'MMMM',
                                                       'bin_1&&g5': 'MMKK'}},
                      'bin_2': {'PF00001.masked.faa': {'bin_2&&g1': 'ACDE',
                                                       'bin_2&&g2': 'ACDE'}},
                      'bin_3': {}}

        self.aaiDir = os.path.join(self.outDir, 'storage', DefaultValues.AAI_QA_DIR)
        for binId, files in alignments.items():
            os.makedirs(os.path.join(self.outDir, 'bins', binId))
            for f, seqs in files.items():
                os.makedirs(os.path.join(self.aaiDir, binId), exist_ok=True)
                writeFasta(seqs, os.path.join(self.aaiDir, binId, f))

        # alignments of a bin that is not in the output directory are ignored
        os.makedirs(os.path.join(self.aaiDir, 'bin_x'))
        writeFasta({'bin_x&&g1': 'AAAA', 'bin_x&&g2': 'AAAA'}, os.path.join(self.aaiDir, 'bin_x', 'PF00001.masked.faa'))

    def tearDown(self):
        shutil.rmtree(self.outDir)

    def runAai(self):
        alignmentFile = os.path.join(self.outDir, 'alignments.txt')
        aai = AminoAcidIdentity()
        aai.run(0.9, self.outDir, alignmentFile)
        with open(alignmentFile) as f:
            alignmentOutput = f.read()

        return aai, alignmentOutput

    def testArchiveAndDirectory(self):
        # CheckM <= 1.2.x output
        aaiDir, alignmentsDir = self.runAai()
        self.assertEqual(set(aaiDir.aaiRawScores), {'bin_1', 'bin_2'})
        self.assertEqual(len(aaiDir.aaiRawScores['bin_1']['PF00001']), 3)
        self.assertAlmostEqual(aaiDir.aaiRawScores['bin_1']['PF00001'][0], 0.9)
        self.assertEqual(aaiDir.aaiRawScores['bin_2']['PF00001'], [1.0])

        # CheckM >= 1.3.0 output
        archiveDirectory(self.aaiDir, os.path.join(self.outDir, 'storage', DefaultValues.AAI_QA_ARCHIVE))
        self.assertFalse(os.path.exists(self.aaiDir))
        aaiArchive, alignmentsArchive = self.runAai()

        self.assertEqual(aaiArchive.aaiRawScores, aaiDir.aaiRawScores)
        self.assertEqual(aaiArchive.aaiHetero, aaiDir.aaiHetero)
        self.assertEqual(aaiArchive.aaiMeanBinHetero, aaiDir.aaiMeanBinHetero)
        self.assertEqual(alignmentsArchive, alignmentsDir)


class VerifyMarkerFileResolution(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(prefix='checkm_test_')
        self.plain = os.path.join(self.tmpDir, 'lineage.ms')
        self.compressed = self.plain + '.gz'

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def _write(self, path, mtime):
        with openFile(path, 'wt') as f:
            f.write('marker file\n')
        os.utime(path, (mtime, mtime))

    def testResolveMarkerFile(self):
        """The old name finds the compressed marker file; if both exist, the newest is used."""
        from checkm.main import OptionsParser
        parser = OptionsParser()

        self._write(self.compressed, 2000)
        self.assertEqual(parser._resolveMarkerFile(self.plain), self.compressed)

        # e.g., a CheckM <= 1.2.x output directory reused with CheckM >= 1.3.0
        self._write(self.plain, 1000)
        with self.assertLogs('timestamp', level='WARNING'):
            self.assertEqual(parser._resolveMarkerFile(self.plain), self.compressed)
            self.assertEqual(parser._resolveMarkerFile(self.compressed), self.compressed)

        os.utime(self.plain, (3000, 3000))
        with self.assertLogs('timestamp', level='WARNING'):
            self.assertEqual(parser._resolveMarkerFile(self.plain), self.plain)

        # custom HMM marker files are used as given
        hmmFile = os.path.join(self.tmpDir, 'custom.hmm')
        self._write(hmmFile, 1000)
        self.assertEqual(parser._resolveMarkerFile(hmmFile), hmmFile)


if __name__ == "__main__":
    unittest.main()
