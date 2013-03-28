#!/usr/bin/env python
###############################################################################
#                                                                             #
#    database.py                                                              #
#                                                                             #
#    Interface to a database of taxon specific markers                        #
#                                                                             #
#    Copyright (C) Connor Skennerton                                          #
#                                                                             #
###############################################################################
#                                                                             #
#                888b     d888  .d8888b.   .d8888b.  888    d8P               #
#                8888b   d8888 d88P  Y88b d88P  Y88b 888   d8P                #
#                88888b.d88888 888    888        888 888  d8P                 #
#                888Y88888P888 888             .d88P 888d88K                  #
#                888 Y888P 888 888         .od888P"  8888888b                 #
#                888  Y8P  888 888    888 d88P"      888  Y88b                #
#                888   "   888 Y88b  d88P 888"       888   Y88b               #
#                888       888  "Y8888P"  888888888  888    Y88b              #
#                                                                             #
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
__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2012, 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.3.2"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

import sys
import sqlite3
import os
import re
from simplehmmer.hmmmodelparser import HmmModelParser, HmmModel
from collections import defaultdict
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Annotation(object):
    """Representation of a Pfam annotation from IMG"""
    def __init__(self, Id, Name=None, Length=None, Count=0):
        super(Annotation, self).__init__()
        self.pfam_id = Id
        self.pfam_name = Name
        self.pfam_length = int(Length)
        self.count = count


class Marker(object):
    def __init__(self):
        super(Marker, self).__init__()
        self.model = None
        self.name = None
        self.dbid = None
        self.taxon = []
        self.taxon_count = []
        self.taxon_single_copy = []

    def taxon_data(self, tax, count, single):
        self.taxon.append(tax)
        self.taxon_count.append(count)
        self.taxon_single_copy.append(single)  


class MarkerDB(object):
    ''' Interface to the database of markers
    '''
    taxonomy_string_re = re.compile('[;,:]?\s?[kpcofgs]__')

    def __init__(self, db=None):
        ''' Set up the database connection
        '''
        self.con = db
        if db is not None:
            self.con = sqlite3.connect(db)
    
    def _parseTaxonomyString(self, taxonomy):
        tax = self.taxonomy_string_re.split(taxonomy)[1:]
        if tax[0] == '':
            return tax[1:]
        return tax

    def _getTaxonIdFromTaxString(self, taxString):
        args = self._parseTaxonomyString(taxString)
        cursor = self.con.cursor()
        ranks = ['Domain', 'Phylum', 'Class', '"Order"', 'Family', 'Genus', 'Species']
        query = []
        for rank, value in zip(ranks, args):
            query.append(' %s = \'%s\' ' % (rank, value))
        query.append(' %s IS NULL' % ranks[len(args)])
        query_string = 'AND'.join(query)
        result = cursor.execute('SELECT Id FROM taxons WHERE %s' % query_string)
        return result.fetchall()

    def _getUniversalMarkers(self, singleCopy=True):
        cur = self.con.cursor()
        if singleCopy:
            return cur.execute('''SELECT bacteria_single_markers.Model FROM bacteria_single_markers
                                    INNER JOIN archaea_single_markers 
                                    ON bacteria_single_markers.Id=archaea_single_markers.Id'''
                                )
        return cur.execute('''SELECT bacteria_markers.Model FROM bacteria_markers
                                    INNER JOIN archaea_markers 
                                    ON bacteria_markers.Id=archaea_markers.Id'''
                                )

    def getModelsForTaxon(self, taxString, singleCopy=True):
        ''' Return a string containing Hmmer3 formatted models

            wrap an SQL query which could be something like:
            SELECT Model FROM markers WHERE Id IN (
                SELECT Marker FROM marker_mapping WHERE Taxon = (
                    SELECT Id FROM taxons WHERE Domain = 'Bacteria' AND Phylum = 'Proteobacteria' AND Class = 'Betaproteobacteria'
                    )
                );
        '''
        if taxString == 'universal':
            return self._getUniversalMarkers(singleCopy=singleCopy)

        tax_id = self._getTaxonIdFromTaxString(taxString)
        if len(tax_id) > 1:
            raise RuntimeError("Taxon string returns more than one lineage "\
                "please be more specific")
        try:
            tax_id = tax_id[0]
        except IndexError, e:
            print tax_id
            raise e

        cur = self.con.cursor()

        if singleCopy:
            result = cur.execute('''SELECT Model FROM markers WHERE Id IN (SELECT Marker FROM
                marker_mapping WHERE SingleCopy = 1 AND Taxon = ?)''', tax_id)
        else:
            result = cur.execute('''SELECT Model FROM markers WHERE Id IN (SELECT Marker FROM
                marker_mapping WHERE Taxon = ?)''', tax_id)

        return result

    def generateModelFiles(self, taxString, outFile, singleCopy=True):
        results = self.getModelsForTaxon(taxString, singleCopy=singleCopy)
        with open(outFile, 'w') as fp:
            for result in results:
                fp.write(result[0])


    #----------------------------
    # Database generation code below

    def makeTables(self):
        self.con.execute("""CREATE TABLE taxons("Id" INTEGER PRIMARY KEY AUTOINCREMENT,
                        "Domain" TEXT, 
                        "Phylum" TEXT, 
                        "Class" TEXT, 
                        "Order" TEXT, 
                        "Family" TEXT, 
                        "Genus" TEXT, 
                        "Species" TEXT,
                        "Count" INTEGER
                        );""")

        self.con.execute("""CREATE TABLE markers("Id" INTEGER PRIMARY KEY AUTOINCREMENT,
                        Name TEXT,
                        Acc TEXT,
                        Model TEXT
                        );""")

        self.con.execute("""CREATE TABLE marker_mapping(Id INTEGER PRIMARY KEY AUTOINCREMENT, 
                                Taxon INTEGER, 
                                Marker INTEGER, 
                                SingleCopy INTEGER,
                                "Count" INTEGER,
                                FOREIGN KEY(Marker) REFERENCES markers(Id), 
                                FOREIGN KEY(Taxon) REFERENCES taxons(Id)
                                );""")
        self.con.execute("""CREATE VIEW bacteria_markers AS SELECT * FROM markers WHERE Id IN (
                                SELECT Marker FROM marker_mapping WHERE Taxon IN (
                                    SELECT Id FROM taxons WHERE (Domain = 'Bacteria') AND Phylum IS NULL
                                )
                            );"""
                        )
        self.con.execute("""CREATE VIEW archaea_markers AS SELECT * FROM markers WHERE Id IN (
                                SELECT Marker FROM marker_mapping WHERE Taxon IN (
                                    SELECT Id FROM taxons WHERE (Domain = 'Archaea') AND Phylum IS NULL
                                )
                            );"""
                        )
        self.con.execute("""CREATE VIEW bacteria_single_markers AS SELECT * FROM markers WHERE Id IN (
                                SELECT Marker FROM marker_mapping WHERE SingleCopy = 1 AND Taxon IN (
                                    SELECT Id FROM taxons WHERE (Domain = 'Bacteria') AND Phylum IS NULL
                                )
                            );"""
                        )
        self.con.execute("""CREATE VIEW archaea_single_markers AS SELECT * FROM markers WHERE Id IN (
                                SELECT Marker FROM marker_mapping WHERE SingleCopy = 1 AND Taxon IN (
                                    SELECT Id FROM taxons WHERE (Domain = 'Archaea') AND Phylum IS NULL
                                )
                            );"""
                        )

    def addMarker(self, cur, model, overwrite=False):
        ''' take a simpleHmmer model object and add it to markers
        '''
        cur.execute('''SELECT Id, Name, Acc FROM markers WHERE Name = ? AND Acc = ?''',
                (model.name, model.acc))
        result = cur.fetchone()
        if result:
            if overwrite:
                cur.execute('''UPDATE markers SET Model=? WHERE Name=? AND Acc=?''',
                        (str(model), model.name, model.acc))
                return result[0]
            else:
                return result[0]
                #raise RuntimeError("A marker with Name=%s and Acc=%s already exists in"\
                #    "the database.  Use -f to force overwriting" % (model.name, model.acc))
        else:
            cur.execute('''INSERT INTO markers (Name, Acc, Model) VALUES (?, ?, ?);''', (model.name, model.acc,
                str(model)))
            return cur.lastrowid

    def addTaxon(self, cur, taxonRanks, genomeCount):
        ''' Add in a taxonomy  to the database
            taxonRanks must be a tuple with 7 elements
        '''
        if not isinstance(taxonRanks, tuple) or len(taxonRanks) != 7:
            raise RuntimeError('The taxon ranks must be a tuple with 7 elements')

        cur.execute('''SELECT Id FROM taxons 
                    WHERE Domain = ? AND
                    Phylum = ? AND
                    Class = ? AND
                    "Order" = ? AND
                    Family = ? AND
                    Genus = ? AND
                    Species = ?''', taxonRanks)
        result = cur.fetchone()

        if not result:
            input_data = list(taxonRanks)
            input_data.append(genomeCount)
            cur.execute('''INSERT INTO taxons (Domain, Phylum, Class, "Order",
                    Family, Genus, Species, "Count") VALUES (?, ?, ?, ?, ?, ?, ?, ?)''',
                    tuple(input_data))
            
            return cur.lastrowid
        else:
            return result[0]

    def addMarkerToTaxon(self, cur, taxonId, markerId, Count, singleCopy=True):
        cur.execute('''INSERT INTO marker_mapping (Taxon, Marker, "Count", SingleCopy)
                VALUES (?, ?, ?, ?)''', (taxonId, markerId, Count, int(singleCopy)))
        return cur.lastrowid

    def _parseMetadata(self, metadataFile):
        '''Take the IMG metadata file and return the taxonomies and counts
        '''
        data = defaultdict(list)
        genomes = []
        with open(metadataFile) as fh:
            for line in fh:
                fields = line.split('\t')
                if fields[2] == 'Finished' and (fields[1] == 'Bacteria' or fields[1] == 'Archaea'):
                    genomes.append(fields[0])
                    taxonomy = []
                    taxonomy.append(fields[1])
                    taxonomy.extend(fields[6:11])
                    max_length = len(taxonomy)+1
                    for i in range(1,max_length):
                        remainder = [None] * (max_length - i)
                        tax = taxonomy[0:i]
                        tax.extend(remainder)
                        tax = tuple(tax)
                        data[tax].append(fields[0])
        return data, genomes


    def _getPfamMarkersFromGenomeAnnotation(self, annotationFile):
        ''' Figure out the pfams that are in this genome
            annotationFile is a python file object
        '''
        counts = {}
        for line in annotationFile:
            fields = line.split('\t')
            if fields[8] == 'pfam_id':
                continue
            try:
                counts[fields[8]] += 1
            except KeyError:
                counts[fields[8]] = 1 
        return counts
            
    def calculateConservedMarkers(self, markerCounts, totalGenomes, minimumConservation, verbose=False):
        conserved_single_markers = {}
        conserved_all_markers = {}
        for marker_name, (single_count, total_count) in markerCounts.items():
            #(single_count, total_count) = data
            single_fraction = float(single_count) / float(totalGenomes)
            total_fraction = float(total_count) / float(totalGenomes)
            #if verbose:
            #    print marker_name, single_fraction, total_fraction, single_count, total_count, totalGenomes
            if single_fraction >= minimumConservation:
                conserved_single_markers['PF' + marker_name[4:]] = single_count
                if verbose:
                    print marker_name, single_fraction, total_fraction

            if total_fraction >= minimumConservation:
                conserved_all_markers['PF' + marker_name[4:]] = total_count
                if verbose:
                    print marker_name, single_fraction, total_fraction

        return conserved_single_markers, conserved_all_markers

    def getMarkersForTaxonomy(self, perGenomeCounts, genomesInTaxonomy, minimumConservation=0.95, verbose=False):
        ''' figure out which markers are common in a taxonomy
            and with are single copy or duplicated
        '''
        markers_in_taxonomy = {}
        if verbose:
            print '\tGenomes:',len(genomesInTaxonomy)
        for genome in genomesInTaxonomy:
            cur_genome_count = perGenomeCounts[genome]
            for marker_name, count in cur_genome_count.items():
                if count == 1:
                    try:
                        markers_in_taxonomy[marker_name][0] += 1
                    except KeyError:
                        markers_in_taxonomy[marker_name] = [1, 0]

                # The count of the marker for this taxonomy regardless of copy number
                try:
                    markers_in_taxonomy[marker_name][1] += 1
                except KeyError:
                    markers_in_taxonomy[marker_name] = [0, 1]
        return self.calculateConservedMarkers(markers_in_taxonomy, len(genomesInTaxonomy), minimumConservation, verbose=verbose)

    def makeDB(self, args):
        genome_taxon_mapping, genomes = self._parseMetadata(args.metadata)

        per_genome_counts = {}
        for genome in genomes:
            annotations = os.path.join(args.img_base, genome, genome + '.pfam.tab.txt')
            if not os.path.exists(annotations):
                print 'Cannot find %s' % annotations
                continue
            per_genome_counts[genome] = self._getPfamMarkersFromGenomeAnnotation(open(annotations))

        if self.con is None:
            self.con = sqlite3.connect(args.database)
        
        cur = self.con.cursor()
        if args.make:
            self.makeTables()
            self.con.commit()

        all_markers = set()
        marker_map = {}
        #sorted_taxons = sorted(genome_taxon_mapping.keys())
        #for t in sorted_taxons:
        #    print t
        for taxon, genomes in genome_taxon_mapping.items():    
            if len(genomes) < args.min_genome:
                continue
            
            verbose = False
            taxon_id = self.addTaxon(cur, taxon, len(genomes))
            self.con.commit()
            # if taxon_id == 151:
            #     verbose = True
            # if verbose:
            #     print taxon
            cons_single, cons_all = self.getMarkersForTaxonomy(per_genome_counts, genomes, args.min_cons, verbose=verbose)
            all_markers |= set(cons_all.keys())
            for marker in cons_all.keys():
                try:
                    marker_map[marker].taxon_data(taxon_id, cons_all[marker], marker in cons_single)
                except KeyError:
                    marker_map[marker] = Marker()
                    marker_map[marker].name = marker
                    marker_map[marker].taxon_data(taxon_id, cons_all[marker], marker in cons_single)
                #all_taxon_markers[taxon][marker] = marker in cons_single

        # extract all models from pfam
        parser = HmmModelParser(args.pfam)
        for model in parser.parse():
            name = model.acc.split('.')
            if name[0] in all_markers:
                marker_db_id = self.addMarker(cur, model)
                marker = marker_map[name[0]]
                for i in range(len(marker.taxon)):
                    self.addMarkerToTaxon(cur, marker.taxon[i], marker_db_id,
                            marker.taxon_count[i], marker.taxon_single_copy[i])

        self.con.commit()
