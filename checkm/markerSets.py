###############################################################################
#
# markerSets.py - Calculate and process marker sets.
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
import re
import sys
import ast
import logging
import functools
import pickle
import gzip

from checkm.defaultValues import DefaultValues
from checkm.hmmer import HMMERRunner
from checkm.hmmerModelParser import HmmModelParser

from checkm.util.pfam import PFAM


# patterns for parsing string representations of sets of marker genes, e.g.
# {'PF00001.1', 'TIGR00001'}, set(['PF00001.1']), set([]) and set()
_SET_STR_RE = re.compile(r"""\{([^{}]*)\}|set\(\[([^\[\]]*)\]\)|set\(\)""")
_QUOTED_STR_RE = re.compile(r"'([^'\\]*)'" + r'|"([^"\\]*)"')


def _parseQuotedItems(content):
    """Parse comma-separated quoted strings, returning None if they are not in this format."""

    # fast path for items without commas, quotes or escapes within them
    items = []
    for token in content.split(','):
        token = token.strip()
        if not token:
            continue

        quote = token[0]
        if len(token) >= 2 and quote in '\'"' and token[-1] == quote and quote not in token[1:-1] and '\\' not in token:
            items.append(token[1:-1])
        else:
            items = None
            break

    if items is not None:
        return items

    # general case
    if _QUOTED_STR_RE.sub('', content).strip(' ,'):
        return None

    return [q.group(1) if q.group(1) is not None else q.group(2) for q in _QUOTED_STR_RE.finditer(content)]


@functools.lru_cache(maxsize=None)
def _parseSetStrs(s):
    """Parse a comma-separated sequence of set representations.

    Returns a tuple of frozensets, or None if the string is not in one of the
    formats written by CheckM. Results are cached since marker sets are repeated
    many times (e.g., the marker sets of each lineage in a lineage marker file).
    """
    sets = []
    pos = 0
    for m in _SET_STR_RE.finditer(s):
        if s[pos:m.start()].strip(' ,'):
            return None

        content = m.group(1) if m.group(1) is not None else (m.group(2) or '')
        items = _parseQuotedItems(content)
        if items is None:
            return None

        sets.append(frozenset(items))
        pos = m.end()

    if s[pos:].strip(' ,'):
        return None

    return tuple(sets)


def parseMarkerSetStr(s):
    """Parse the string representation of a marker set (a list of sets of marker genes).

    This replaces the use of eval(), which was slow and evaluated file contents as code.
    The sets are returned as frozensets which may be shared between marker sets, since
    marker sets are never modified in place (they are replaced when refined).
    """
    s = s.strip()
    sets = None
    if s.startswith('[') and s.endswith(']'):
        sets = _parseSetStrs(s[1:-1])

    if sets is None:
        # fall back to Python's (safe) literal parser for any other formatting
        return [frozenset(x) for x in ast.literal_eval(s)]

    return list(sets)


def parseSetStr(s):
    """Parse the string representation of a single set of marker genes."""
    s = s.strip()
    sets = _parseSetStrs(s)
    if sets is None or len(sets) != 1:
        return set(ast.literal_eval(s))

    return set(sets[0])


class BinMarkerSets():
    """A collection of one or more marker sets associated with a bin."""

    # type of marker set
    TAXONOMIC_MARKER_SET = 1
    TREE_MARKER_SET = 2
    HMM_MODELS_SET = 3

    def __init__(self, binId, markerSetType):
        self.logger = logging.getLogger('timestamp')
        self.markerSets = []
        self.binId = binId
        self.markerSetType = markerSetType
        self.selectedLinageSpecificMarkerSet = None

    def addMarkerSet(self, markerSet):
        """Add marker set to bin."""
        self.markerSets.append(markerSet)

    def markerSetIter(self):
        """Generator function for iterating over marker sets."""
        for markerSet in self.markerSets:
            yield markerSet

    def getMarkerGenes(self):
        """Get marker genes from all marker sets."""
        markerGenes = set()
        for ms in self.markerSets:
            markerGenes.update(ms.getMarkerGenes())

        return markerGenes

    def mostSpecificMarkerSet(self):
        return self.markerSets[0]

    def selectedMarkerSet(self):
        """Return the 'selected' marker set for this bin."""
        if self.markerSetType == self.TAXONOMIC_MARKER_SET:
            return self.mostSpecificMarkerSet()
        elif self.markerSetType == self.TREE_MARKER_SET:
            return self.selectedLinageSpecificMarkerSet
        else:
            # there should be a single marker set associated with this bin
            if len(self.markerSets) == 1:
                return self.markerSets[0]

            self.logger.error('Expected a single marker set to be associated with each bin.')
            sys.exit(1)

    def setLineageSpecificSelectedMarkerSet(self, selectedMarkerSetMap):
        """Select the lineage-specific marker set used to assess the bin.

        selectedMarkerSetMap maps the UID of each node in the genome tree to
        the UID of the marker set selected for genomes placed below that node.
        The marker sets of the bin are ordered from most specific to the root.
        """
        uidToMarkerSet = {}
        for ms in self.markerSets:
            uidToMarkerSet.setdefault(ms.UID, ms)

        self.selectedLinageSpecificMarkerSet = None
        for ms in self.markerSets:
            # This is a hack for the reduced tree. Since not all
            # marker sets are in the reduced tree it is possible the
            # selected marker set might not be available. In this case,
            # follow the selected marker sets until one associated with
            # the bin is found. If the chain ends (e.g., at a marker set
            # that selects itself) without finding one, move to the next
            # (less specific) marker set of the bin. Ideally, this could be
            # avoided by just forcing in the selected marker set.
            selectedId = selectedMarkerSetMap.get(ms.UID)
            visited = set()
            while selectedId is not None and selectedId not in visited:
                if selectedId in uidToMarkerSet:
                    self.selectedLinageSpecificMarkerSet = uidToMarkerSet[selectedId]
                    return

                visited.add(selectedId)
                selectedId = selectedMarkerSetMap.get(selectedId)

        # something has gone wrong
        self.logger.error('Failed to set a selected lineage-specific marker set for bin %s.' % self.binId)
        sys.exit(1)

    def removeMarkers(self, markersToRemove):
        """Remove specified markers from all marker sets."""
        for markerSet in self.markerSets:
            curMarkerSet = markerSet.markerSet

            newMarkerSet = []
            for ms in curMarkerSet:
                newMS = ms - markersToRemove

                if len(newMS) != 0:
                    newMarkerSet.append(newMS)

            markerSet.markerSet = newMarkerSet

    def write(self, fout):
        """Write marker set to file."""
        fout.write(self.binId)
        fout.write('\t' + str(len(self.markerSets)))
        for ms in self.markerSets:
            fout.write('\t' + str(ms))
        fout.write('\n')

    def read(self, line):
        """Construct bin marker set data from line."""
        lineSplit = line.split('\t')
        numMarkerSets = int(lineSplit[1])
        for i in range(0, numMarkerSets):
            uid = lineSplit[i * 4 + 2]
            lineageStr = lineSplit[i * 4 + 3]
            numGenomes = int(lineSplit[i * 4 + 4])
            markerSet = parseMarkerSetStr(lineSplit[i * 4 + 5])
            self.markerSets.append(MarkerSet(uid, lineageStr, numGenomes, markerSet))


class MarkerSet():
    """A collection of marker genes organized into co-located sets."""
    def __init__(self, UID, lineageStr, numGenomes, markerSet):
        self.logger = logging.getLogger('timestamp')

        self.UID = UID  # unique ID of marker set
        self.lineageStr = lineageStr  # taxonomic string associated with marker set
        self.numGenomes = numGenomes  # number of genomes used to calculate marker set
        self.markerSet = markerSet  # marker genes organized into co-located sets

    def __repr__(self):
        # marker genes are written in sorted order so marker files are reproducible
        # (the string must remain a valid Python list of sets; see parseMarkerSetStr)
        markerSetStr = '[' + ', '.join([self._setStr(m) for m in self.markerSet]) + ']'
        return str(self.UID) + '\t' + self.lineageStr + '\t' + str(self.numGenomes) + '\t' + markerSetStr

    @staticmethod
    def _setStr(markers):
        """String representation of a set with its elements in sorted order."""
        if len(markers) == 0:
            return 'set()'

        return '{' + ', '.join([repr(m) for m in sorted(markers)]) + '}'

    def size(self):
        """Number of marker genes and marker gene sets."""
        numMarkerGenes = 0
        for m in self.markerSet:
            numMarkerGenes += len(m)

        return numMarkerGenes, len(self.markerSet)

    def numMarkers(self):
        """Number of marker genes."""
        return self.size()[0]

    def numSets(self):
        """Number of marker sets."""
        return len(self.markerSet)

    def getMarkerGenes(self):
        """Get marker genes within marker set."""
        markerGenes = set()
        for m in self.markerSet:
            for marker in m:
                markerGenes.add(marker)

        return markerGenes
        
    def removeMarkers(self, markersToRemove):
        """Remove specified markers from marker sets."""
        newMarkerSet = []
        for ms in self.markerSet:
            newMS = ms - markersToRemove

            if len(newMS) != 0:
                newMarkerSet.append(newMS)

        self.markerSet = newMarkerSet

    def genomeCheck(self, hits, bIndividualMarkers):
        """Calculate genome completeness and contamination."""
        if self.numMarkers() == 0:
            # e.g., all markers were removed with --exclude_markers
            if not getattr(self, 'bWarnedNoMarkers', False):
                self.logger.warning('Marker set %s (%s) contains no marker genes; reporting 0%% completeness and contamination.' % (str(self.UID), self.lineageStr))
                self.bWarnedNoMarkers = True
            return 0.0, 0.0

        if bIndividualMarkers:
            present = 0
            multiCopyCount = 0
            for marker in self.getMarkerGenes():
                count = len(hits.get(marker, []))
                if count > 0:
                    present += 1
                    multiCopyCount += (count - 1)

            percComp = 100 * float(present) / self.numMarkers()
            percCont = 100 * float(multiCopyCount) / self.numMarkers()
        else:
            comp = 0.0
            cont = 0.0
            numSets = 0
            for ms in self.markerSet:
                if len(ms) == 0:
                    continue

                numSets += 1
                present = 0
                multiCopy = 0
                for marker in ms:
                    count = len(hits.get(marker, []))
                    if count == 1:
                        present += 1
                    elif count > 1:
                        present += 1
                        multiCopy += (count - 1)

                comp += float(present) / len(ms)
                cont += float(multiCopy) / len(ms)

            percComp = 100 * comp / numSets
            percCont = 100 * cont / numSets

        return percComp, percCont


class MarkerSetParser():
    """Parse marker set file."""

    def __init__(self, threads=1):
        self.logger = logging.getLogger('timestamp')
        self.numThreads = threads

    def getMarkerSets(self, outDir, binIds, markerFile, excludeMarkersFile=None):
        """Determine marker set for each bin."""

        # determine type of marker set file
        markerFileType = self.markerFileType(markerFile)

        # get marker set for each bin
        binIdToBinMarkerSets = {}

        if markerFileType == BinMarkerSets.TAXONOMIC_MARKER_SET:
            binMarkerSets = self.parseTaxonomicMarkerSetFile(markerFile)

            for binId in binIds:
                binIdToBinMarkerSets[binId] = binMarkerSets
        elif markerFileType == BinMarkerSets.TREE_MARKER_SET:
            binIdToBinMarkerSets = self.parseLineageMarkerSetFile(markerFile)
        else:
            markers = [set()]
            modelParser = HmmModelParser(markerFile)
            for model in modelParser.parse():
                markers[0].add(model.acc)
            markerSet = MarkerSet(0, "N/A", -1, markers)

            for binId in binIds:
                binMarkerSets = BinMarkerSets(binId, BinMarkerSets.HMM_MODELS_SET)
                binMarkerSets.addMarkerSet(markerSet)
                binIdToBinMarkerSets[binId] = binMarkerSets

        # remove marker genes specified by user or marked for exclusion
        markersToExclude = set()
        if excludeMarkersFile:
            markersToExclude = self.readExcludeMarkersFile(excludeMarkersFile)

        markersToExclude.update(DefaultValues.MARKERS_TO_EXCLUDE)
        for binId, binMarkerSet in binIdToBinMarkerSets.items():
            binMarkerSet.removeMarkers(markersToExclude)

        return binIdToBinMarkerSets

    def readExcludeMarkersFile(self, excludeMarkersFile):
        """Parse file specifying markers to exclude."""
        markersToExclude = set()
        for line in open(excludeMarkersFile):
            if line[0] == '#':
                continue

            marker = line.strip()
            markersToExclude.add(marker)

        return markersToExclude

    def markerFileType(self, markerFile):
        """Determine type of marker file."""
        with open(markerFile, 'r') as f:
            header = f.readline()

        if DefaultValues.TAXON_MARKER_FILE_HEADER in header:
            return BinMarkerSets.TAXONOMIC_MARKER_SET
        elif DefaultValues.LINEAGE_MARKER_FILE_HEADER in header:
            return BinMarkerSets.TREE_MARKER_SET
        elif 'HMMER3' in header:
            return BinMarkerSets.HMM_MODELS_SET
        else:
            self.logger.error('Unrecognized file type: ' + markerFile)
            sys.exit(1)

    def hmmAccessions(self, binMarkerSet, pfam):
        """Get accessions of the HMMs to search for a marker set.

        These are the marker genes in the marker set along with all genes
        from the same Pfam clan as any marker gene.
        """
        markerGenes = binMarkerSet.getMarkerGenes()
        return markerGenes | pfam.genesInSameClan(markerGenes)

    def fetchHmms(self, accessions, outputFile):
        """Create a HMM file containing the specified models from the CheckM HMM file."""

        # accessions are sorted so the order of models in the HMM file
        # (and hence downstream output) is deterministic
        keyFile = outputFile + '.keys'
        with open(keyFile, 'w') as fout:
            for modelAcc in sorted(accessions):
                fout.write(modelAcc + '\n')

        try:
            HMMERRunner(mode='fetch').fetch(DefaultValues.HMM_MODELS, keyFile, outputFile, bKeyFile=True)
        finally:
            os.remove(keyFile)

    def parseTaxonomicMarkerSetFile(self, markerSetFile):
        """Parse marker set from a taxonomic-specific marker set file."""
        with open(markerSetFile) as f:
            f.readline()  # skip header

            binLine = f.readline()
            taxonId = binLine.split('\t')[0]
            binMarkerSets = BinMarkerSets(taxonId, BinMarkerSets.TAXONOMIC_MARKER_SET)
            binMarkerSets.read(binLine)

        return binMarkerSets

    def parseLineageMarkerSetFile(self, markerSetFile):
        """Parse marker sets from a lineage-specific marker set file."""

        # read all marker sets
        selectedMarkerSetMap = self.parseSelectedMarkerSetMap()

        binIdToBinMarkerSets = {}
        with open(markerSetFile) as f:
            f.readline()  # skip header

            for line in f:
                lineSplit = line.split('\t')
                binId = lineSplit[0]

                binMarkerSets = BinMarkerSets(binId, BinMarkerSets.TREE_MARKER_SET)
                binMarkerSets.read(line)

                # determine selected marker set
                binMarkerSets.setLineageSpecificSelectedMarkerSet(selectedMarkerSetMap)

                binIdToBinMarkerSets[binId] = binMarkerSets

        return binIdToBinMarkerSets

    def parseSelectedMarkerSetMap(self):
        selectedMarkerSetMap = {}
        for line in open(DefaultValues.SELECTED_MARKER_SETS):
            lineSplit = line.split('\t')
            internalID = lineSplit[0]
            selectedID = lineSplit[1].rstrip()

            selectedMarkerSetMap[internalID] = selectedID

        return selectedMarkerSetMap

    def writeBinModels(self, binIdToModels, filename):
        """Save HMM model info for each bin to file."""

        self.logger.info('Saving HMM info to file.')

        with gzip.open(filename, 'wb') as output:
            pickle.dump(binIdToModels, output, pickle.HIGHEST_PROTOCOL)

    def loadBinModels(self, filename):
        """Read HMM model info for each bin from file."""

        self.logger.info('Reading HMM info from file.')

        with gzip.open(filename, 'rb') as f:
            binIdToModels = pickle.load(f)

        return binIdToModels
