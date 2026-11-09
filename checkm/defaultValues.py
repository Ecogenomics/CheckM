###############################################################################
#
# defaultValues.py - store default values used in many places in CheckM
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
from checkm.checkmData import DBManager


class DefaultValues():
    """Default values for filenames and common constants."""

    __DBM = DBManager()

    # set of markers recognized to be unreliable. These are often
    # ubiquitous, single-copy genes, but ones which are challenging
    # to correctly annotate with the PFAM and TIGRFAM models.
    MARKERS_TO_EXCLUDE = {'TIGR00398', 'TIGR00399'}

    E_VAL = 1e-10
    LENGTH = 0.7
    PSEUDOGENE_LENGTH = 0.3

    TAXON_MARKER_FILE_HEADER = '# [Taxon Marker File]'
    LINEAGE_MARKER_FILE_HEADER = '# [Lineage Marker File]'

    SEQ_CONCAT_CHAR = '&&'

    CHECKM_DATA_DIR = __DBM.config.values["dataRoot"]
    PHYLO_HMM_MODELS = os.path.join(
        CHECKM_DATA_DIR, 'hmms', 'phylo.hmm')
    HMM_MODELS = os.path.join(CHECKM_DATA_DIR, 'hmms', 'checkm.hmm')
    PFAM_CLAN_FILE = os.path.join(CHECKM_DATA_DIR, 'pfam', 'Pfam-A.hmm.dat')

    SELECTED_MARKER_SETS = os.path.join(
        CHECKM_DATA_DIR, 'selected_marker_sets.tsv')
    TAXON_MARKER_SETS = os.path.join(CHECKM_DATA_DIR, 'taxon_marker_sets.tsv')

    GENOME_TREE_DIR = os.path.join(CHECKM_DATA_DIR, 'genome_tree')
    # pplacer reference packages (used by gPlace)
    PPLACER_REF_PACKAGE_FULL = os.path.join(
        GENOME_TREE_DIR, 'genome_tree_full.refpkg')
    GENOME_TREE = 'genome_tree.tre'
    GENOME_TREE_FASTA = 'genome_tree.fasta'
    GENOME_TREE_DEREP = 'genome_tree.derep.txt'
    GENOME_TREE_TAXONOMY = 'genome_tree.taxonomy.tsv'
    GENOME_TREE_METADATA = 'genome_tree.metadata.tsv'
    GENOME_TREE_MISSING_DUPLICATE = 'missing_duplicate_genes_50.tsv'
    DISTRIBUTION_DIR = os.path.join(CHECKM_DATA_DIR, 'distributions')

    PHYLO_HMM_MODEL_INFO = 'phylo_hmm_info.pkl.gz'
    CHECKM_HMM_MODEL_INFO = 'checkm_hmm_info.pkl.gz'

    # output files are gzip-compressed from CheckM 1.3.0; files are read with
    # checkm.util.compression.openFile(), which also reads the uncompressed
    # files (e.g., genes.faa rather than genes.faa.gz) of earlier versions
    HMMER_TABLE_PHYLO_OUT = 'hmmer.tree.txt.gz'
    HMMER_PHYLO_OUT = 'hmmer.tree.ali.txt.gz'

    HMMER_TABLE_OUT = 'hmmer.analyze.txt.gz'
    HMMER_OUT = 'hmmer.analyze.ali.txt.gz'

    PRODIGAL_AA = 'genes.faa.gz'
    PRODIGAL_NT = 'genes.fna.gz'
    PRODIGAL_GFF = 'genes.gff.gz'

    # placement files (gPlace replaced pplacer in CheckM 1.3.0; file names are unchanged for compatibility, apart from compression)
    PPLACER_CONCAT_SEQ_OUT = 'concatenated.fasta.gz'
    PPLACER_JSON_OUT = 'concatenated.pplacer.json.gz'
    PPLACER_OUT = 'pplacer.out'
    PPLACER_TREE_OUT = 'concatenated.tre.gz'

    BIN_STATS_PHYLO_OUT = 'bin_stats.tree.tsv.gz'
    BIN_STATS_OUT = 'bin_stats.analyze.tsv.gz'

    BIN_STATS_EXT_OUT = 'bin_stats_ext.tsv.gz'
    MARKER_GENE_STATS = 'marker_gene_stats.tsv.gz'

    # alignments of markers with multiple hits in a bin (written by analyze, read by qa);
    # CheckM <= 1.2.x wrote these to the directory storage/aai_qa
    AAI_QA_DIR = 'aai_qa'
    AAI_QA_ARCHIVE = 'aai_qa.tar.gz'

    CONTIG_BREAK = 'NNNNNNNNNN'

    UNBINNED = 'unbinned'

    # quality report written by lineage_wf and taxonomy_wf (output of 'qa -o 2 --tab_table')
    QUALITY_REPORT = 'quality_report.tsv'

    MIN_SEQ_LEN_GC_STD = 1000

