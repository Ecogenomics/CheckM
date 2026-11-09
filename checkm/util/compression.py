###############################################################################
#
# compression.py - read and write gzip-compressed output files
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

"""Read and write gzip-compressed files.

CheckM >= 1.3.0 writes most output files gzip-compressed (e.g., genes.faa.gz),
while earlier versions wrote them uncompressed (e.g., genes.faa). Files are
read with resolvePath()/openFile(), which accept either name, so output
directories of earlier versions can still be used.

Files are compressed with xopen, which uses ISA-L (python-isal) when available.
Compression level 1 is used as it is fast and compresses these files almost
as well as higher levels. Files are compressed in the calling process
(threads=0) as CheckM already processes bins in parallel.
"""

import os
import shutil
import tarfile

from xopen import xopen

COMPRESS_LEVEL = 1


def compressedName(path):
    """Name of a path with a '.gz' extension (added if not already present)."""
    if path.endswith('.gz'):
        return path

    return path + '.gz'


def uncompressedName(path):
    """Name of a path without a '.gz' extension."""
    if path.endswith('.gz'):
        return path[0:-3]

    return path


def resolvePath(path):
    """Path of an existing file, allowing for compression.

    Returns path if it exists, otherwise path + '.gz' if it exists,
    otherwise the path without its '.gz' extension if it exists (e.g., an
    uncompressed output file of CheckM <= 1.2.x). If none of these exist,
    path is returned unchanged.
    """
    if os.path.exists(path):
        return path

    if os.path.exists(path + '.gz'):
        return path + '.gz'

    if path.endswith('.gz') and os.path.exists(path[0:-3]):
        return path[0:-3]

    return path


def openFile(path, mode='rt'):
    """Open a file that may be gzip-compressed.

    Files opened for reading are located with resolvePath() and are
    decompressed if they are gzip-compressed (determined from their
    content, so uncompressed files can be read regardless of their name).
    Files opened for writing are gzip-compressed if path ends with '.gz'.
    """
    if 'r' in mode:
        path = resolvePath(path)
        if isGzipped(path):
            return xopen(path, mode, threads=0, format='gz')

        return open(path, mode)

    return xopen(path, mode, compresslevel=COMPRESS_LEVEL, threads=0)


def isGzipped(path):
    """Check if a file is gzip-compressed (from its content rather than its name)."""
    with open(path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def compressFile(src, dst=None, remove=True):
    """Compress a file (src may itself be gzip-compressed).

    The file is written to dst (default: src + '.gz') and src is removed
    if remove is True. Returns the path of the compressed file.
    """
    if dst is None:
        dst = src + '.gz'

    with openFile(src, 'rb') as fin, xopen(dst, 'wb', compresslevel=COMPRESS_LEVEL, threads=0, format='gz') as fout:
        shutil.copyfileobj(fin, fout, 1024 * 1024)

    if remove:
        os.remove(src)

    return dst


def hasContent(path):
    """Check if a file, which may be gzip-compressed, exists and contains data.

    The size of a gzip-compressed file can't be used for this, as a
    compressed file without any data is ~20 bytes.
    """
    path = resolvePath(path)
    if not os.path.exists(path):
        return False

    with openFile(path, 'rb') as f:
        return len(f.read(1)) != 0


def archiveDirectory(inputDir, archiveFile, remove=True):
    """Create a gzip-compressed tar archive of the files in a directory.

    Members are named by their path relative to inputDir and are added in
    sorted order. An existing archive is replaced. The directory is
    removed if remove is True.
    """
    tmpArchiveFile = archiveFile + '.tmp'
    with tarfile.open(tmpArchiveFile, 'w:gz', compresslevel=COMPRESS_LEVEL) as tar:
        for root, dirs, files in os.walk(inputDir):
            dirs.sort()
            for f in sorted(files):
                filePath = os.path.join(root, f)
                tar.add(filePath, arcname=os.path.relpath(filePath, inputDir), recursive=False)
    os.replace(tmpArchiveFile, archiveFile)

    if remove:
        shutil.rmtree(inputDir)
