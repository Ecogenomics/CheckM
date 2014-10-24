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

	__DBM = DBManager()

	"""Default values for filenames and common constants."""
	E_VAL = 1e-10
	LENGTH = 0.7

	TAXON_MARKER_FILE_HEADER = '# [Taxon Marker File]'
	LINEAGE_MARKER_FILE_HEADER = '# [Lineage Marker File]'

	SEQ_CONCAT_CHAR = '&&'

	CHECKM_DATA_DIR = __DBM.config.values["dataRoot"]
	PHYLO_HMM_MODELS = phyloHMMs = os.path.join(CHECKM_DATA_DIR, 'hmms', 'phylo.hmm')
	HMM_MODELS = os.path.join(CHECKM_DATA_DIR, 'hmms', 'checkm.hmm')
	PFAM_CLAN_FILE = os.path.join(CHECKM_DATA_DIR, 'pfam', 'Pfam-A.hmm.dat')

	IMG_METADATA_FILE = os.path.join(CHECKM_DATA_DIR,'img', 'img_metadata.tsv')
	REDUNDANT_TIGRFAM_FILE = os.path.join(CHECKM_DATA_DIR,'pfam', 'tigrfam2pfam.tsv')

	SELECTED_MARKER_SETS = os.path.join(CHECKM_DATA_DIR, 'selected_marker_sets.tsv')
	TAXON_MARKER_SETS = os.path.join(CHECKM_DATA_DIR, 'taxon_marker_sets.tsv')

	GENOME_TREE_DIR = os.path.join(CHECKM_DATA_DIR, 'genome_tree')
	PPLACER_REF_PACKAGE = os.path.join(GENOME_TREE_DIR, 'genome_tree_prok.refpkg')
	GENOME_TREE_FINAL = 'genome_tree.final.tre'
	DISTRIBUTION_DIR = os.path.join(CHECKM_DATA_DIR, 'distributions')

	PHYLO_HMM_MODEL_INFO = 'phylo_hmm_info.pkl.gz'
	CHECKM_HMM_MODEL_INFO = 'checkm_hmm_info.pkl.gz'

	HMMER_TABLE_PHYLO_OUT = 'hmmer.tree.txt'
	HMMER_PHYLO_OUT = 'hmmer.tree.ali.txt'

	HMMER_TABLE_OUT = 'hmmer.analyze.txt'
	HMMER_OUT = 'hmmer.analyze.ali.txt'

	PRODIGAL_AA = 'prodigal.faa'
	PRODIGAL_NT = 'prodigal.fna'
	PRODIGAL_GFF = 'prodigal.gff'

	PPLACER_CONCAT_SEQ_OUT = 'concatenated.fasta'
	PPLACER_JSON_OUT = 'concatenated.pplacer.json'
	PPLACER_OUT = 'pplacer.out'
	PPLACER_TREE_OUT = 'concatenated.tre'

	BIN_STATS_PHYLO_OUT = 'bin_stats.tree.tsv'
	SEQ_STATS_PHYLO_OUT = 'seq_stats.tree.tsv'

	BIN_STATS_OUT = 'bin_stats.analyze.tsv'
	SEQ_STATS_OUT = 'seq_stats.analyze.tsv'

	BIN_STATS_EXT_OUT = 'bin_stats_ext.tsv'
	MARKER_GENE_STATS = 'marker_gene_stats.tsv'

	CONTIG_BREAK = 'NNNNNNNNNN'

	UNBINNED = 'unbinned'

	MIN_SEQ_LEN_GC_STD = 1000
