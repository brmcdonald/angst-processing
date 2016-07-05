import sys
import os
import re
import math
from scipy import stats
from Bio import Phylo
import numpy as np
sys.path.append('/home/bradon/scripts/BRM_LIB')
import brm_stdlib as brm
import brm_seqlib as seq

#
#Environmental Variables for KEGG
#
k2g = "/home/bradon/dbs/KEGG/ko2gene.tbl"
kegg_pathway = "/home/bradon/dbs/KEGG/Pathway_KOs.list"

#Kegg categories that are worth analyzing in bacterial genomes
GOOD_CAT = {"Metabolism_of_Terpenoids_and_Polyketides":"","Biosynthesis_of_Other_Secondary_Metabolites":"",
"Cell_Growth_and_Death":"","Amino_Acid_Metabolism":"","Carbohydrate_Metabolism":"","Energy_Metabolism":"",
"Glycan_Biosynthesis_and_Metabolism":"","Lipid_Metabolism":"","Membrane_Transport":"",
"Metabolism_of_Cofactors_and_Vitamins":"","Metabolism_of_Other_Amino_Acids":"","Nucleotide_Metabolism":"",
"Replication_and_Repair":"","Signal_Transduction":"","Transcription":"","Translation":"",
"Xenobiotics_Biodegradation_and_Metabolism":"", "Mobile_Elements":""}

debug = 0

if debug == 1:
	GOOD_CAT["Test1"] = ""
	GOOD_CAT["Test2"] = ""

#Kegg pathways that are worth analyzing in bacterial genomes
GOOD_PATH = {"Basal_transcription_factors":"","RNA_polymerase":"","translation:":"","Aminoacyl-tRNA_biosynthesis":"",
"RNA_transport":"","Ribosome":"","tRNAs":"","signal_trans:":"","Calcium_signaling_pathway":"",
"ErbB_signaling_pathway":"","Jak-STAT_signaling_pathway":"","Phosphatidylinositol_signaling_system":"",
"Plant_hormone_signal_transduction":"","Two-component_system":"","motility:":"","Bacterial_chemotaxis":"",
"Flagellar_assembly":"","vitamins:":"","Biotin_metabolism":"","Cyanoamino_acid_metabolism":"",
"D-Alanine_metabolism":"","D-Arginine_and_D-ornithine_metabolism":"",
"D-Glutamine_and_D-glutamate_metabolism":"","Folate_biosynthesis":"","Glutathione_metabolism":"",
"Lipoic_acid_metabolism":"","Nicotinate_and_nicotinamide_metabolism":"","One_carbon_pool_by_folate":"",
"Pantothenate_and_CoA_biosynthesis":"","Phosphonate_and_phosphinate_metabolism":"",
"Porphyrin_and_chlorophyll_metabolism":"","Retinol_metabolism":"","Riboflavin_metabolism":"",
"Selenocompound_metabolism":"","Taurine_and_hypotaurine_metabolism":"","Thiamine_metabolism":"",
"Ubiquinone_and_other_terpenoid-quinone_biosynthesis":"","Vitamin_B6_metabolism":"",
"beta-Alanine_metabolism":"","glycans:":"","Glycosaminoglycan_biosynthesis_-_chondroitin_sulfate":"",
"Glycosaminoglycan_biosynthesis_-_heparan_sulfate":"",
"Glycosaminoglycan_biosynthesis_-_keratan_sulfate":"","Glycosaminoglycan_degradation":"",
"Glycosphingolipid_biosynthesis_-_ganglio_series":"","Glycosphingolipid_biosynthesis_-_globo_series":"",
"Glycosphingolipid_biosynthesis_-_lacto_and_neolacto_series":"",
"Glycosylphosphatidylinositol(GPI)-anchor_biosynthesis":"","Lipopolysaccharide_biosynthesis":"",
"Mucin_type_O-glycan_biosynthesis":"","N-Glycan_biosynthesis":"","Other_glycan_degradation":"",
"Other_types_of_O-glycan_biosynthesis":"","Peptidoglycan_biosynthesis":"",
"Various_types_of_N-glycan_biosynthesis":"","rep_repair:":"","Base_excision_repair":"","DNA_replication":"",
"Homologous_recombination":"","Mismatch_repair":"","Non-homologous_end-joining":"",
"Nucleotide_excision_repair":"","2ndary_cats:":"","Acridone_alkaloid_biosynthesis":"",
"Anthocyanin_biosynthesis":"","Benzoxazinoid_biosynthesis":"","Betalain_biosynthesis":"",
"Biosynthesis_of_12-,_14-_and_16-membered_macrolides":"","Biosynthesis_of_ansamycins":"",
"Biosynthesis_of_siderophore_group_nonribosomal_peptides":"",
"Biosynthesis_of_type_II_polyketide_backbone":"","Biosynthesis_of_type_II_polyketide_products":"",
"Biosynthesis_of_vancomycin_group_antibiotics":"","Brassinosteroid_biosynthesis":"",
"Butirosin_and_neomycin_biosynthesis":"","Caffeine_metabolism":"","Carotenoid_biosynthesis":"",
"Clavulanic_acid_biosynthesis":"","Diterpenoid_biosynthesis":"","Flavone_and_flavonol_biosynthesis":"",
"Flavonoid_biosynthesis":"","Glucosinolate_biosynthesis":"","Indole_alkaloid_biosynthesis":"",
"Isoflavonoid_biosynthesis":"","Isoquinoline_alkaloid_biosynthesis":"","Monoterpenoid_biosynthesis":"",
"Novobiocin_biosynthesis":"","Penicillin_and_cephalosporin_biosynthesis":"","Phenylpropanoid_biosynthesis":""
,"Polyketide_sugar_unit_biosynthesis":"","Puromycin_biosynthesis":"","Sesquiterpenoid_biosynthesis":"",
"Stilbenoid,_diarylheptanoid_and_gingerol_biosynthesis":"","Streptomycin_biosynthesis":"",
"Terpenoid_backbone_biosynthesis":"","Tetracycline_biosynthesis":"",
"Tropane,_piperidine_and_pyridine_alkaloid_biosynthesis":"","Zeatin_biosynthesis":"",
"beta-Lactam_resistance":"","lipids:":"","Arachidonic_acid_metabolism":"",
"Biosynthesis_of_unsaturated_fatty_acids":"","Ether_lipid_metabolism":"",
"Fatty_acid_biosynthesis":"","Fatty_acid_elongation_in_mitochondria":"",
"Fatty_acid_metabolism":"","Glycerolipid_metabolism":"","Glycerophospholipid_metabolism":"",
"Linoleic_acid_metabolism":"","Primary_bile_acid_biosynthesis":"","Secondary_bile_acid_biosynthesis":"",
"Sphingolipid_metabolism":"","Steroid_biosynthesis":"","Steroid_hormone_biosynthesis":"",
"Synthesis_and_degradation_of_ketone_bodies":"","alpha-Linolenic_acid_metabolism":"",
"nucleotide_metab:":"","Purine_metabolism":"","Pyrimidine_metabolism":"","cell_growth:":"",
"Apoptosis":"","Cell_cycle":"","Cell_cycle_-_Caulobacter":"","carb_metab:":"",
"Amino_sugar_and_nucleotide_sugar_metabolism":"","Ascorbate_and_aldarate_metabolism":"",
"Butanoate_metabolism":"","C5-Branched_dibasic_acid_metabolism":"","Citrate_cycle_(TCA_cycle)":"",
"Fructose_and_mannose_metabolism":"","Galactose_metabolism":"","Glycolysis_/_Gluconeogenesis":"",
"Glyoxylate_and_dicarboxylate_metabolism":"","Inositol_phosphate_metabolism":"",
"Pentose_and_glucuronate_interconversions":"","Pentose_phosphate_pathway":"","Propanoate_metabolism":"",
"Pyruvate_metabolism":"","Starch_and_sucrose_metabolism":"","membrane_trans:":"","ABC_transporters":"",
"Bacterial_secretion_system":"","Phosphotransferase_system_(PTS)":"",
"xenobiotics:":"","1,1,1-Trichloro-2,2-bis(4-chlorophenyl)ethane_(DDT)_degradation":"",
"Aminobenzoate_degradation":"","Atrazine_degradation":"","Benzoate_degradation":"",
"Bisphenol_degradation":"","Caprolactam_degradation":"","Chloroalkane_and_chloroalkene_degradation":"",
"Chlorocyclohexane_and_chlorobenzene_degradation":"","Dioxin_degradation":"",
"Drug_metabolism_-_cytochrome_P450":"","Drug_metabolism_-_other_enzymes":"","Ethylbenzene_degradation":"",
"Fluorobenzoate_degradation":"","Metabolism_of_xenobiotics_by_cytochrome_P450":"",
"Naphthalene_degradation":"","Nitrotoluene_degradation":"","Polycyclic_aromatic_hydrocarbon_degradation":"",
"Styrene_degradation":"","Toluene_degradation":"","Xylene_degradation":"","aa_metab:":"",
"Alanine,_aspartate_and_glutamate_metabolism":"","Arginine_and_proline_metabolism":"",
"Cysteine_and_methionine_metabolism":"","Glycine,_serine_and_threonine_metabolism":"",
"Histidine_metabolism":"","Lysine_biosynthesis":"","Lysine_degradation":"","Phenylalanine_metabolism":"",
"Phenylalanine,_tyrosine_and_tryptophan_biosynthesis":"","Tryptophan_metabolism":"","Tyrosine_metabolism":"",
"Valine,_leucine_and_isoleucine_biosynthesis":"","Valine,_leucine_and_isoleucine_degradation":"",
"energy_metab:":"","Carbon_fixation_in_photosynthetic_organisms":"",
"Carbon_fixation_pathways_in_prokaryotes":"","Methane_metabolism":"","Nitrogen_metabolism":"",
"Oxidative_phosphorylation":"","Photosynthesis":"","Photosynthesis_-_antenna_proteins":"",
"Sulfur_metabolism":""}

class Path():
	def __init__(self, name):
		self.name = name
		self.ko = {}		#Dict of KOs in the pathway/category

		#Dicts of event data: {NODE:event_number}
		self.lgt_source = {}
		self.lgt_target = {}
		self.lgt = {}
		self.dup = {}
		self.los = {}
		self.spc = {}
		self.brn = {}
		self.total_gene = {}
		self.total_event = {}

		#Total events for the category across the tree
		self.lgt_global = 0.0
		self.dup_global = 0.0
		self.los_global = 0.0
		self.spc_global = 0.0
		self.brn_global = 0.0
		self.total_gene_global = 0.0
		self.total_event_global = 0.0

#Contains data for each node of the species tree, simpler than full version in data parser
class Node():
	def __init__(self, obj, flag):
		self.term_list = {}		#List of terminal nodes ie genomes
		self.dir_dec = {}		#List of direct decendent nodes
		self.dec = {}			#List of all decendent nodes,includes self
		self.node_obj = obj		#Node object from biopython
		self.parent = ''		#Parent node
		self.genome = flag		#Whether or not the node is a genome (terminal)
		self.los = {}			#Number of loss events occured at node for a cat
		self.dup = {}			#Number of duplication events occured at node for a cat
		self.brn = {}			#Genes born at node
		self.spc = {}			#Genes derived from strict vertical inheritence
		self.event = {}			#Number of events involving a cat

		self.los_total = 0		#Total loss events
		self.dup_total = 0		#Total duplication events
		self.brn_total = 0		#Total genes born at node
		self.spc_total = 0		#Total speciation events
		self.event_total = 0	#Total events
		self.angst_name = ''	#Angst name, based on terminal node names
		self.name = ''			#Name assigned by auto namer

		self.total = 0.0			#Total number of genes at a node
		self.cats = {}				#Number of genes in each cat

		self.lgt_source = {}		#LGT events from current node to [key]
		self.lgt_target = {}		#LGT events to current node from [key]
		self.lgt_total_source = 0.0	#Total lgt events from current node
		self.lgt_total_target = 0.0	#Total lgt events to current node
		self.lgt_source_gene = {}	#Number of lgt events involving each gene family
		self.lgt_target_gene = {}	#Number of lgt events involving each gene family

		self.include = {}			#Only used for combining data from multiple nodes

class Gene():
	def __init__(self):
		self.name = ''
		self.brn = ''
		self.genome = ''
		self.fam = ''
		self.dup = []
		self.lgt = []	#Each event is a list: [source,destination]
		self.spc = []


####################################################################################################
#######################################    Functions     ###########################################
####################################################################################################


def import_data(LEAF_parse,NODE_N_parse,FAMS_parse,tree_parse):
	"""Used to move global variables from main parsing script"""
	global LEAF
	global NODE_N
	global FAMS
	global tree
	LEAF = LEAF_parse
	NODE_N = NODE_N_parse
	FAMS = FAMS_parse
	tree = tree_parse
	return

def r_fet(a,b,c,d,name):
	"""Returns text required to run FET in R"""

	fet_input = ("""
	over <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
	pval <- fisher.test(over, alternative=\"greater\")$p.value
	odds <- fisher.test(over, alternative=\"greater\")$estimate
	fet_data_over <- rbind(fet_data_over, c("{4}",pval,odds))

	under <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
	pval <- fisher.test(under, alternative=\"less\")$p.value
	odds <- fisher.test(under, alternative=\"less\")$estimate
	fet_data_under <- rbind(fet_data_under, c("{4}",pval,odds))

	both <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
	pval <- fisher.test(both, alternative=\"two.sided\")$p.value
	odds <- fisher.test(both, alternative=\"two.sided\")$estimate
	fet_data_both <- rbind(fet_data_both, c("{4}",pval,odds))\n""".format(a,b,c,d,name))

	return fet_input

def piechart_tree_proportion(NODE_N,outfile,fam='ALL'):
	"""Generates piechart for iTOL showing the proportion of each event type"""
	out = open(outfile, 'w')
	most_events = 0.0
	OUT_NODE = {}

	if fam == "ALL":
		out.write("LABELS\tLGT\tDUP\tLOS\tBRN\tSPC\nCOLORS\t#E00000\t#3DB319\t#2460C7\t#FFFF05\t#858282\n")
		TOTAL = {}
		for i in NODE_N:
			total = 0.0
			total += NODE_N[i].dup_total
			total += NODE_N[i].lgt_total_target
			total += NODE_N[i].los_total
			total += NODE_N[i].brn_total
			total += NODE_N[i].spc_total
			TOTAL[NODE_N[i].name] = total

			if total > most_events:
				most_events = total

			if total > 0:
				dup = (float(NODE_N[i].dup_total) / total) * 100
				lgt = (float(NODE_N[i].lgt_total_target) / total) * 100
				los = (float(NODE_N[i].los_total) / total) * 100
				brn = (float(NODE_N[i].brn_total) / total) * 100
				spc = (float(NODE_N[i].spc_total) / total) * 100
				OUT_NODE[NODE_N[i].name] = "{0}\tR100\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(NODE_N[i].name,lgt,dup,los,brn,spc)
		#Need to calculate piechart size by total events at node / most events at a node
		for i in OUT_NODE:
			out.write(OUT_NODE[i])

	else:
		out.write("LABELS\tLGT\tDUP\tLOS\tBRN\tSPC\nCOLORS\t#E00000\t#3DB319\t#2460C7\t#FFFF05\t#858282\n")
		TOTAL = {}
		for i in NODE_N:
			total = 0.0
			total += NODE_N[i].dup[fam]
			total += NODE_N[i].lgt_target_gene[fam]
			total += NODE_N[i].los[fam]
			if NODE_N[i].brn.has_key(fam):
				total += 1
			total += NODE_N[i].spc[fam]
			TOTAL[NODE_N[i].name] = total

			if total > most_events:
				most_events = total

			if total > 0:
				dup = (float(NODE_N[i].dup[fam]) / total) * 100
				lgt = (float(NODE_N[i].lgt_target_gene[fam]) / total) * 100
				los = (float(NODE_N[i].los[fam]) / total) * 100
				spc = (float(NODE_N[i].spc[fam]) / total) * 100
				if NODE_N[i].brn.has_key(fam):
					brn = (float(NODE_N[i].brn[fam]) / total) * 100
				else:
					brn = 0

				OUT_NODE[NODE_N[i].name] = "{0}\tR100\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(NODE_N[i].name,lgt,dup,los,brn,spc)
		#Need to calculate piechart size by total events at node / most events at a node
		for i in OUT_NODE:
			out.write(OUT_NODE[i])

def lgt_tree(outfile):
	"""Generates LGT for iTOL"""
	out = open(outfile, 'w')
	for i in NODE_N:
		for target in NODE_N[i].lgt_source:
			out.write("{0}\t{1}\t#2460C7\t{2}\n"
			.format(NODE_N[i].name,NODE_N[target].name,target))
	out.close()

def cat_data_table(outfile,build_flag):
	"""Generates a table of events for each gene family"""
	global NODE_N
	global FAMS

	out = open(outfile, 'w')
	out.write("Family\tName\tLGT\tDUP\tLOSS\tSPC\tFinal_Count\n")
	if build_flag == "N":
		for fam in FAMS:
			out.write("{6}\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
			.format(KEGG_NAME[fam], FAMS[fam].lgt, FAMS[fam].dup, FAMS[fam].los,FAMS[fam].spc,FAMS[fam].total,fam))
		out.close()
	else:
		DATA = build_kegg("k2g")
		KEGG_NAME = DATA[2]

		for fam in FAMS:
			out.write("{6}\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
			.format(KEGG_NAME[fam], FAMS[fam].lgt, FAMS[fam].dup, FAMS[fam].los,FAMS[fam].spc,FAMS[fam].total,fam))
		out.close()


def node_isolation(NODE_N,out_file,out_type,border_node='',tree=''):
	"""Determines the proportion of LGT events that occur from outside a clade"""
	PRINT_DATA = {}
	RET_DATA = {}
	if out_type == 'itol':
		out_sum = open("{0}_summary.out".format(out_file),'w')
		out_tbl = open("{0}_tbl.out".format(out_file),'w')
		out_itol = open('{0}.itol'.format(out_file),'w')

	#Recording 
	if border_node:
		dec = NODE_N[border_node].dec
		border_dec = {}
		border_dec[NODE_N[border_node].name] = ""
		for j in dec:
			border_dec[NODE_N[j].name] = ""

	for i in NODE_N:
		PRINT_DATA[i] = {}
		RET_DATA[i] = {}
		#Count LGT from subclade, out of subclade, external, and from each other node
		PRINT_DATA[i]["out_lgt"] = 0.0
		PRINT_DATA[i]["in_clade_lgt"] = 0.0
		PRINT_DATA[i]["in_border_lgt"] = 0.0
		PRINT_DATA[i]["total_lgt"] = 0.0
		PRINT_DATA[i]["external"] = 0.0
		PRINT_DATA[i]["vertical"] = 0.0

		RET_DATA[i]['in_clade_lgt'] = {}
		RET_DATA[i]['in_border_lgt'] = {}
		RET_DATA[i]['out_lgt'] = {}
		RET_DATA[i]['brn'] = ''

		for j in NODE_N:
			PRINT_DATA[i][j] = 0.0

	for i in NODE_N:
		read = []		#List of nodes included in this clade
		read.append(i)
		dec = NODE_N[i].dec
		dec_angst = {}
		dec_angst[NODE_N[i].name] = ""
		for j in dec:
			read.append(j)
			dec_angst[NODE_N[j].name] = ""

		#For each node in the subclade
		for j in read:

			#For each LGT event into the node
			for lgt in NODE_N[j].lgt_target:
				PRINT_DATA[i]["total_lgt"] += len(NODE_N[j].lgt_target[lgt])

				if dec_angst.has_key(lgt):
					PRINT_DATA[i]["in_clade_lgt"] += len(NODE_N[j].lgt_target[lgt])
					RET_DATA[i]['in_clade_lgt'][lgt] = NODE_N[j].lgt_target[lgt]
				elif border_dec.has_key(lgt) and border_node:
					PRINT_DATA[i]['in_border_lgt'] += len(NODE_N[j].lgt_target[lgt])
					RET_DATA[i]['in_border_lgt'][lgt] = NODE_N[j].lgt_target[lgt]
				else:
					PRINT_DATA[i]["out_lgt"] += len(NODE_N[j].lgt_target[lgt])
					RET_DATA[i]['out_lgt'][lgt] = NODE_N[j].lgt_target[lgt]
					if PRINT_DATA[i][NODE_N[lgt].name] == 0:
						PRINT_DATA[i][NODE_N[lgt].name] += len(NODE_N[j].lgt_target[lgt])

			#Assume gene birth events are LGT from outside the tree
			for brn in NODE_N[j].brn:
				RET_DATA[i]['brn'] = NODE_N[j].brn
				PRINT_DATA[i]["out_lgt"] += 1
				PRINT_DATA[i]["total_lgt"] += 1
			PRINT_DATA[i]["vertical"] += NODE_N[j].spc_total

	if out_type == 'itol':
		out_sum.write("Node\tin_lgt\tout_lgt\ttotal\tperc_in\n")
		for i in PRINT_DATA:
			if NODE_N[i].genome == 0:
				perc = PRINT_DATA[i]["in_clade_lgt"] / PRINT_DATA[i]["total_lgt"]
				out_sum.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(i,PRINT_DATA[i]["in_clade_lgt"],PRINT_DATA[i]["out_lgt"],PRINT_DATA[i]["total_lgt"], perc))
		out_sum.close()

		out_tbl.write("NODE")
		for i in sorted(PRINT_DATA):
			out_tbl.write("\t{0}".format(i))
		out_tbl.write("\n")
		for i in sorted(PRINT_DATA):
			out_tbl.write(i)
			for j in sorted(PRINT_DATA[i]):
				if PRINT_DATA.has_key(j):
					out_tbl.write("\t{0}".format(PRINT_DATA[i][j]))
			out_tbl.write("\n")

		out_itol.write("LABELS\tLGT_IN\tLGT_Border\tLGT_OUT\tVERTICAL\nCOLORS\t#34A134\t#257FE6\t#E62F12\t#A1A1A1\n")
		for i in PRINT_DATA:
			if NODE_N[i].genome == 0:
				total_genes = PRINT_DATA[i]["total_lgt"] + PRINT_DATA[i]["vertical"]
				perc_vertical = PRINT_DATA[i]["vertical"] / total_genes
				perc_in = PRINT_DATA[i]["in_clade_lgt"] / total_genes
				perc_border = PRINT_DATA[i]['in_border_lgt'] / total_genes
				perc_external = PRINT_DATA[i]["out_lgt"] / total_genes
				out_itol.write("{0}\tR100\t{1}\t{2}\t{3}\t{4}\n".format(i,perc_in,perc_border,perc_external,perc_vertical))
		out_sum.close()
		out_itol.close()

	if tree:
		out = open('test','w')
		from Bio import Phylo
		ISO_DIST = {}
		for node in NODE_N['NODE92'].dec:
			if NODE_N[node].genome == 0:
				perc_iso = (PRINT_DATA[node]["in_clade_lgt"] / PRINT_DATA[node]["total_lgt"])

				max_dist = 0
				for i in NODE_N[node].dec:
					if NODE_N[i].genome == 1:
						for j in NODE_N[node].dec:
							if NODE_N[j].genome == 1:
								if tree.distance(i,j) > max_dist:
									max_dist = (tree.distance(i,j))
				out.write('{0}\t{1}\t{2}\n'.format(node,perc_iso,max_dist))
		out.close()

	if out_type == 'itol':
		return [RET_DATA,PRINT_DATA]

def build_kegg(flag):
	"""Generates lookup tables for KEGG data"""
	global NODE_N
	global FAMS
	LOOKUP = {}
	LOOKUP_GENE = {}

	#Kegg gene name lookup table: {KO:gene_name}
	if flag == "k2g":
		kegg = open(k2g, 'r')
		kegg = kegg.readlines()
		for ln in kegg:
			ln = ln.replace("\n","")
			ln = ln.split("\t")
			LOOKUP_GENE[ln[0]] = ln[1]
		for fam in FAMS:
			if not LOOKUP_GENE.has_key(fam):
				LOOKUP_GENE[fam] = "{0} family protein".format(fam)

	#Kegg category assignments: {KO:[categories]}
	elif flag == "k2c":
		kegg = open (kegg_pathway, 'r')
		kegg = kegg.readlines()
		for ln in kegg:
			ln = ln.replace("\n","")
			if "|" in ln:
				ln = ln.split(": ")
				if len(ln) == 2:
					ln[1] = ln[1].replace("|","")
					cat = ln[1].replace(" ","_")
				else:
					cat = ''
				if not GOOD_CAT.has_key(cat):
					cat = ''
			elif cat != '':
				ln = ln.split(".")
				for i in ln:
					if i != '':
						if LOOKUP.has_key(i):
							LOOKUP[i][cat] = ""
						else:
							LOOKUP[i] = {}
							LOOKUP[i][cat] = ""

	#Kegg pathway assignments: {KO:[pathways]}
	elif flag == "k2p":
		kegg = open (kegg_pathway, 'r')
		kegg = kegg.readlines()
		for ln in kegg:
			ln = ln.replace("\n","")
			if "|" in ln:
				ln = ln.split(": ")
				if len(ln) == 3:
					ln[2] = ln[2].replace("|","")
					cat = ln[2].replace(" ","_")
				else:
					cat = ''
				if not GOOD_PATH.has_key(cat):
					cat = ''
			elif cat != '':
				ln = ln.split(".")
				for i in ln:
					if i != '':
						if LOOKUP.has_key(i):
							LOOKUP[i][cat] = ""
						else:
							LOOKUP[i] = {}
							LOOKUP[i][cat] = ""

	else:
		print "Cannot understand tag {0}\nexiting".format(flag)
		sys.exit()

	CATS = {}
	if flag == "k2c" or flag == "k2p":
		cat_names = []
		for cat in LOOKUP:
			for i in LOOKUP[cat]:
				if cat:
					cat_names.append(i)

		for cat in cat_names:
			cat_obj = Path(cat)
			CATS[cat] = cat_obj

		for cat in CATS:

			#Make a list of gene fams in the category
			for i in LOOKUP:
				if LOOKUP[i].has_key(cat):
					CATS[cat].ko[i] = ""

			#Go through each node and add the number of events affecting genes in the category
			for node in NODE_N:
				name = NODE_N[node].name
				CATS[cat].dup[name] = 0.0
				CATS[cat].lgt[name] = 0.0
				CATS[cat].los[name] = 0.0
				CATS[cat].spc[name] = 0.0
				CATS[cat].brn[name] = 0.0
				CATS[cat].total_gene[name] = 0.0
				CATS[cat].total_event[name] = 0.0

				#Gene totals are counted from every node in the tree, not just terminal nodes
				for fam in NODE_N[node].cats:
					if CATS[cat].ko.has_key(fam):
						CATS[cat].total_gene[name] += NODE_N[node].cats[fam]
						CATS[cat].total_gene_global += NODE_N[node].cats[fam]

				for fam in NODE_N[node].dup:
					if CATS[cat].ko.has_key(fam):
						CATS[cat].dup[name] += NODE_N[node].dup[fam]
						CATS[cat].total_event[name] += NODE_N[node].dup[fam]
				for fam in NODE_N[node].los:
					if CATS[cat].ko.has_key(fam):
						CATS[cat].los[name] += NODE_N[node].los[fam]
						CATS[cat].total_event[name] += NODE_N[node].los[fam]
				for fam in NODE_N[node].spc:
					if CATS[cat].ko.has_key(fam):
						CATS[cat].spc[name] += NODE_N[node].spc[fam]
						CATS[cat].total_event[name] += NODE_N[node].spc[fam]
				for fam in NODE_N[node].lgt_target_gene:
					if CATS[cat].ko.has_key(fam):
						CATS[cat].lgt[name] += NODE_N[node].lgt_target_gene[fam]
						CATS[cat].total_event[name] += NODE_N[node].lgt_target_gene[fam]
				for fam in NODE_N[node].brn:
					if CATS[cat].ko.has_key(fam):
						CATS[cat].brn[name] += NODE_N[node].brn[fam]
						CATS[cat].total_event[name] += NODE_N[node].brn[fam]

		for i in FAMS:
			if LOOKUP.has_key(i):
				for cat in LOOKUP[i]:
					CATS[cat].lgt_global += FAMS[i].lgt
					CATS[cat].spc_global += FAMS[i].spc
					CATS[cat].dup_global += FAMS[i].dup
					CATS[cat].los_global += FAMS[i].los
					CATS[cat].total_event_global += FAMS[i].event

	return_data = []
	return_data.append(LOOKUP)
	return_data.append(CATS)
	return_data.append(LOOKUP_GENE)
	return return_data

#Run FET
def fet_cat(CATS,total_lgt,total_gene_treewide,event):
	global NODE_N
	#Idiot proofing
	if event == 'hgt':
		event = 'lgt'

	fet_input = open("fet_input_{0}.R".format(event), 'w')

	#Set up output matrices
	fet_input.write("""fet_data_over <- matrix(c(\"Node\",\"FET_P\", \"FET_OR\"), ncol=3)
fet_data_under <- matrix(c(\"Node\",\"FET_P\", \"FET_OR\"), ncol=3)
fet_data_both <- matrix(c(\"Node\",\"FET_P\", \"FET_OR\"), ncol=3)""")

	R_name_lst = "c(\"FET\""
	for cat in CATS:
		R_name_lst += ",\"{0}\"".format(cat)
		a_all = 0
		b_all = 0
		c_all = 0
		d_all = 0

		#For each node
		for node in NODE_N:
			R_name_lst += ",\"{0}\"".format(cat)
			node_name = NODE_N[node].name

			#Is an event more/less likely to occur in a particular node, for each category

			#	Matrix layout:
			#	FOR GENES IN A CATEGORY:
			#
			#	X axis: event
			#	Y axis: node
			#
			#	a	c
			#	b	d
			if event == "lgt":
				a = CATS[cat].lgt[node_name]
				b = CATS[cat].total_gene[node_name] - CATS[cat].lgt[node_name]
				c = CATS[cat].lgt_global - CATS[cat].lgt[node_name]
				d = CATS[cat].total_gene_global - CATS[cat].lgt_global - b
			elif event == "dup":
				a = CATS[cat].dup[node_name]
				b = CATS[cat].total_gene[node_name] - CATS[cat].dup[node_name]
				c = CATS[cat].dup_global - CATS[cat].dup[node_name]
				d = CATS[cat].total_gene_global - CATS[cat].dup_global - b
			elif event == "los":
				a = CATS[cat].los[node_name]
				b = CATS[cat].total_gene[node_name] - CATS[cat].los[node_name]
				if b < 0:
					b = 0
				c = CATS[cat].los_global - CATS[cat].los[node_name]
				d = CATS[cat].total_gene_global - CATS[cat].los_global - b
			elif event == "spc":
				a = CATS[cat].spc[node_name]
				b = CATS[cat].total_gene[node_name] - CATS[cat].spc[node_name]
				c = CATS[cat].spc_global - CATS[cat].spc[node_name]
				d = CATS[cat].total_gene_global - CATS[cat].spc_global - b


			#R commands for running FET on the node
			fet_input.write("""
over <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(over, alternative=\"greater\")$p.value
odds <- fisher.test(over, alternative=\"greater\")$estimate
fet_data_over <- rbind(fet_data_over, c("{4}",pval,odds))

under <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(under, alternative=\"less\")$p.value
odds <- fisher.test(under, alternative=\"less\")$estimate
fet_data_under <- rbind(fet_data_under, c("{4}",pval,odds))

both <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(both, alternative=\"two.sided\")$p.value
odds <- fisher.test(both, alternative=\"two.sided\")$estimate
fet_data_both <- rbind(fet_data_both, c("{4}",pval,odds))\n""".format(a,b,c,d,node_name))

		#Runs FET for all nodes combined: Is a category more likely to be <event> 
		#across the tree
		if event == "lgt":
			a_all = CATS[cat].lgt_global
			b_all = CATS[cat].total_gene_global - CATS[cat].lgt_global
			c_all = total_lgt - CATS[cat].lgt_global
			d_all = total_gene_treewide - total_lgt - b_all
		elif event == 'dup':
			a_all = CATS[cat].dup_global
			b_all = CATS[cat].total_gene_global - CATS[cat].dup_global
			c_all = total_dup - CATS[cat].dup_global
			d_all = total_gene_treewide - total_dup - b_all
		elif event == 'los':
			a_all = CATS[cat].los_global
			b_all = CATS[cat].total_gene_global - CATS[cat].los_global
			c_all = total_los - CATS[cat].los_global
			d_all = total_gene_treewide - total_los - b_all
		elif event == 'spc':
			a_all = CATS[cat].spc_global
			b_all = CATS[cat].total_gene_global - CATS[cat].spc_global
			c_all = total_spc - CATS[cat].spc_global
			d_all = total_gene_treewide - total_spc - b_all

		#R commands for running FET on the category across all nodes
		fet_input.write("""\n#Cat is {5}\n
over <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(over, alternative=\"greater\")$p.value
odds <- fisher.test(over, alternative=\"greater\")$estimate
fet_data_over <- rbind(fet_data_over, c("{4}",pval,odds))

under <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(under, alternative=\"less\")$p.value
odds <- fisher.test(under, alternative=\"less\")$estimate
fet_data_under <- rbind(fet_data_under, c("{4}",pval,odds))

both <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(both, alternative=\"two.sided\")$p.value
odds <- fisher.test(both, alternative=\"two.sided\")$estimate
fet_data_both <- rbind(fet_data_both, c("{4}",pval,odds))\n""".format(a_all,b_all,c_all,d_all,"All_Nodes",cat))


	R_name_lst += ")"

	fet_input.write("""
row.names(fet_data_under) <- {0}
row.names(fet_data_over) <- {0}
row.names(fet_data_both) <- {0}
write.table(fet_data_over, \"{1}_fet_over.tbl\", quote=FALSE, row.names=TRUE, col.names=FALSE, sep=\"\\t\")
write.table(fet_data_both, \"{1}_fet_both.tbl\", quote=FALSE, row.names=TRUE, col.names=FALSE, sep=\"\\t\")
write.table(fet_data_under, \"{1}_fet_under.tbl\", quote=FALSE, row.names=TRUE, col.names=FALSE, sep=\"\\t\")""".format(R_name_lst,event))

	fet_input.close()
	os.system("R CMD BATCH fet_input_{0}.R".format(event))
	#os.system("rm -f fet_input_temp.*")


def cat_total_ttest(nodes1,nodes2,CATS,event):
	"""Determine if the total number of LGT events from a particulary category is significantly different in two sets of nodes"""
	from scipy import stats

	SET1 = {}
	SET2 = {}
	if event == 'lgt':
		for i in CATS.keys():
			SET1[i] = []
			SET2[i] = []
		for c in CATS.keys():
			for i in nodes1:
				if CATS[c].lgt.has_key(i):
					SET1[c].append(CATS[c].lgt[i])
				else:
					SET1[c].append(0)
			for i in nodes2:
				if CATS[c].lgt.has_key(i):
					SET2[c].append(CATS[c].lgt[i])
				else:
					SET2[c].append(0)
	elif event == 'dup':
		for i in CATS.keys():
			SET1[i] = []
			SET2[i] = []
		for c in CATS.keys():
			for i in nodes1:
				SET1[c].append(CATS[c].dup[i])
			for i in nodes2:
				SET2[c].append(CATS[c].dup[i])
	if event == 'los':
		for i in CATS.keys():
			SET1[i] = []
			SET2[i] = []
		for c in CATS.keys():
			for i in nodes1:
				SET1[c].append(CATS[c].los[i])
			for i in nodes2:
				SET2[c].append(CATS[c].los[i])

	print "Category\tT_stat\tPval\tMean_set1\tMean_set2\n"
	for cat in SET1:
		t,pval = stats.ttest_ind(SET1[cat],SET2[cat],equal_var=False)
		print "{0}\t{1}\t{2}\t{3}\t{4}".format(cat,t,pval,np.mean(SET1[cat]),np.mean(SET2[cat]))

def cat_proportion_ttest(nodes1,nodes2,CATS,event):
	"""Determine if the proportion of LGT events from a particulary category is significantly different in two sets of nodes"""
	from scipy import stats

	SET1 = {}
	SET2 = {}
	totalLGT1 = []
	totalLGT2 = []

	for i in nodes1:
		totalLGT1.append(0.0)
	for i in nodes2:
		totalLGT2.append(0.0)


	if event == 'lgt':
		for i in CATS.keys():
			SET1[i] = []
			SET2[i] = []
		for c in CATS.keys():
			for idx, i in enumerate(nodes1):
				if CATS[c].lgt.has_key(i):
					SET1[c].append(CATS[c].lgt[i])
					totalLGT1[idx] += CATS[c].lgt[i]
				else:
					SET1[c].append(0)
			for idx, i in enumerate(nodes2):
				if CATS[c].lgt.has_key(i):
					SET2[c].append(CATS[c].lgt[i])
					totalLGT2[idx] += CATS[c].lgt[i]
				else:
					SET2[c].append(0)

	#avoid division by 0 errors
	for idx, i in enumerate(totalLGT1):
		if i == 0:
			totalLGT1[idx] = 1.0
	for idx, i in enumerate(totalLGT2):
		if i == 0:
			totalLGT2[idx] = 1.0

	print "Category\tT_stat\tPval\tMean_set1\tMean_set2"
	for cat in SET1:
		prop1 = [SET1[cat][idx] / totalLGT1[idx] for idx, i in enumerate(nodes1)]
		prop2 = [SET2[cat][idx] / totalLGT2[idx] for idx, i in enumerate(nodes2)]

		t,pval = stats.ttest_ind(prop1,prop2,equal_var=False)
		print "{0}\t{1}\t{2}\t{3}\t{4}".format(cat,t,pval,np.mean(prop1),np.mean(prop2))


def fet_fam(event):
	global NODE_N
	fet_input = open("fet_input_temp.R", 'w')

	#Set up output matrices
	fet_input.write("""fet_data_over <- matrix(c(\"Node\",\"FET_P\", \"FET_OR\"), ncol=3)
fet_data_under <- matrix(c(\"Node\",\"FET_P\", \"FET_OR\"), ncol=3)
fet_data_both <- matrix(c(\"Node\",\"FET_P\", \"FET_OR\"), ncol=3)""")

	R_name_lst = "c(\"FET\""
	for fam in FAMS:
		R_name_lst += ",\"{0}\"".format(fam)
		a_all = 0
		b_all = 0
		c_all = 0
		d_all = 0

		#For each node
		for node in NODE_N:
			R_name_lst += ",\"{0}\"".format(fam)
			node_name = NODE_N[node].name

			#	Matrix layout:
			#	a	c
			#	b	d
			if event == "lgt":
				a = NODE_N[node].lgt_gene[fam]
				b = NODE_N[node].lgt_total - a
				c = NODE_N[node].event[fam] - a
				d = NODE_N[node].event_total - a
			elif event == "dup":
				a = NODE_N[node].dup[fam]
				b = NODE_N[node].dup_total - a
				c = NODE_N[node].event[fam] - a
				d = NODE_N[node].event_total - a
			elif event == "los":
				a = NODE_N[node].los[fam]
				b = NODE_N[node].los_total - a
				c = NODE_N[node].event[fam] - a
				d = NODE_N[node].event_total - a
			elif event == "spc":
				a = NODE_N[node].spc[fam]
				b = NODE_N[node].spc_total - a
				c = NODE_N[node].event_total - a
				d = NODE_N[node].event_total - a

			#Runs FET for all nodes combined
			a_all += a
			b_all += b
			c_all += c
			d_all += d

			#R commands for running FET on the node
			fet_input.write("""
over <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(over, alternative=\"greater\")$p.value
odds <- fisher.test(over, alternative=\"greater\")$estimate
fet_data_over <- rbind(fet_data_over, c("{4}",pval,odds))

under <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(under, alternative=\"less\")$p.value
odds <- fisher.test(under, alternative=\"less\")$estimate
fet_data_under <- rbind(fet_data_under, c("{4}",pval,odds))

both <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(both, alternative=\"two.sided\")$p.value
odds <- fisher.test(both, alternative=\"two.sided\")$estimate
fet_data_both <- rbind(fet_data_both, c("{4}",pval,odds))\n""".format(a,b,c,d,node_name))

		#R commands for running FET on the category across all nodes
		fet_input.write("""
over <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(over, alternative=\"greater\")$p.value
odds <- fisher.test(over, alternative=\"greater\")$estimate
fet_data_over <- rbind(fet_data_over, c("{4}",pval,odds))

under <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(under, alternative=\"less\")$p.value
odds <- fisher.test(under, alternative=\"less\")$estimate
fet_data_under <- rbind(fet_data_under, c("{4}",pval,odds))

both <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
pval <- fisher.test(both, alternative=\"two.sided\")$p.value
odds <- fisher.test(both, alternative=\"two.sided\")$estimate
fet_data_both <- rbind(fet_data_both, c("{4}",pval,odds))\n""".format(a_all,b_all,c_all,d_all,"All_Nodes"))


	R_name_lst += ")"

	fet_input.write("""
row.names(fet_data_under) <- {0}
row.names(fet_data_over) <- {0}
row.names(fet_data_both) <- {0}
write.table(fet_data_over, \"{1}_fet_over.tbl\", quote=FALSE, row.names=TRUE, col.names=FALSE, sep=\"\\t\")
write.table(fet_data_both, \"{1}_fet_both.tbl\", quote=FALSE, row.names=TRUE, col.names=FALSE, sep=\"\\t\")
write.table(fet_data_under, \"{1}_fet_under.tbl\", quote=FALSE, row.names=TRUE, col.names=FALSE, sep=\"\\t\")""".format(R_name_lst,event))

	fet_input.close()
	os.system("R CMD BATCH fet_input_temp.R")

def event_rate(event,out_format,conv=1.0,includeNode=''):
	"""Event per branch length"""
	global NODE_N

	class RateNode():
		def __init__(self):
			self.event_num = 0.0
			self.age = 0.0
			self.rate = 0.0

	from scipy.stats.mstats import mquantiles
	from bisect import bisect_left

	#Will keep all nodes decendant from the 'includeNode'
	if includeNode:
		keep = NODE_N[includeNode].term_list
	DATA = {}
	data_array = []
	colors = ['#0007DB','#194FFF','#4589FF','#69C8FF','#00BD71','#00BD00','#FFAB19','#D40000']

	for i in NODE_N:
		if event == 'lgt':
			event_total = NODE_N[i].lgt_total_target
		elif event == 'dup':
			event_total = NODE_N[i].dup_total
		elif event == 'los':
			event_total = NODE_N[i].los_total
		elif event == 'brn':
			event_total = NODE_N[i].brn_total

		bl = NODE_N[i].node_obj.branch_length

		if NODE_N[i].node_obj.branch_length > 0:
			read = 1
			genomes = NODE_N[i].node_obj.get_terminals()
			for g in genomes:
				if not keep.has_key(g.name):
					read = 0
			if read == 1:
				rateNode = RateNode()
				rateNode.name = NODE_N[i].name
				rateNode.event_num = event_total
				rateNode.age = bl * conv
				rateNode.rate = (event_total / rateNode.age)
				DATA[rateNode.name] = rateNode
				data_array.append(DATA[rateNode.name].rate)

	bounds = mquantiles(data_array, prob=[0.125, 0.25, 0.375,.5,.625,.75,.875,1], alphap=1, betap=1)

	if out_format == 'itol':
		out = open('{0}_rate.itol'.format(event),'w')
		out.write('NODE_ID\tTYPE\tCOLOR\tLABEL\n')
		for i in DATA:
			color = colors[bisect_left(bounds,DATA[i].rate)]
			out.write('{0}\tclade\t{1}\tno_label\n'.format(i,color))
		out.close()
	elif out_format == 'tbl':
		out = open('{0}_rate.tbl'.format(event),'w')
		out.write('Node\tage\tevent_num\tevent_rate\n')
		for i in DATA:
			printData = (DATA[i].name,DATA[i].age,DATA[i].event_num,DATA[i].rate)
			out.write('{0}\t{1}\t{2}\t{3}\n'.format(*printData))
		out.close()


def tip_isolation(clade_division,out_type,out_file):
	global NODE_N
	try:
		clade_data = open(clade_division,'r')
	except IOError:
		return 'Cant read clade division file'

	if out_type == 'raw':
		raw = 1
	elif out_type == 'perc':
		raw = 0
	else:
		return 'Output type is either raw (numbers) or perc (percentages'

	out = open(out_file,'w')
	out.write('Labels\tIn\tOut\tActino\n')
	out.write('COLORS\t#969696\t#219FFF\t#D60909\n')

	clade = clade_data.readlines()
	clade_data.close()

	#Reading clade file
	CLADE = {}
	for ln in clade:
		ln = ln.replace('\n','')
		if '>' in ln:
			clade = ln.replace('>','')
		elif ln:
			g = ln.split('\t')[0]
			if CLADE.has_key(g):
				print 'Clade divisions must be mutually exclusive: {0} {1} {2}'.format(g,clade,CLADE[g])
				sys.exit()
			CLADE[g] = clade

	for n in NODE_N:
		if NODE_N[n].genome != 5:
			if CLADE.has_key(n):
				home_clade = CLADE[n]
				in_clade = 0.0
				out_clade = 0.0
				actino = 0.0
				all_lgt = 0.0

				for lgt in NODE_N[n].lgt_target:
					if NODE_N[lgt].genome == 1:
						all_lgt += len(NODE_N[n].lgt_target[lgt])
						if CLADE.has_key(lgt):
							if CLADE[lgt] == home_clade:
								in_clade += len(NODE_N[n].lgt_target[lgt])
							elif CLADE[lgt] == 'actino':
								actino += len(NODE_N[n].lgt_target[lgt])
							else:
								out_clade += len(NODE_N[n].lgt_target[lgt])
						else:
							out_clade += len(NODE_N[n].lgt_target[lgt])

				if raw == 0:
					p_in = (in_clade / all_lgt) * 100
					p_out = (out_clade / all_lgt) * 100
					p_actino = (actino / all_lgt) * 100
				elif raw == 1:
					p_in = in_clade
					p_out = out_clade
					p_actino = actino
				out.write('{0}\t{1}\t{2}\t{3}\n'.format(n,p_in,p_out,p_actino))
	out.close()


def get_full_sublcade_data(node_name,ancestral_flag):
	"""Combine data from all nodes in a subclade into one node object, or only ancestral nodes"""
	global NODE_N
	out_node = Node(NODE_N[node_name].node_obj,0)
	out_node.name = node_name
	out_node.term_list = NODE_N[node_name].term_list
	out_node.dir_dec = NODE_N[node_name].dir_dec
	out_node.dec = NODE_N[node_name].dec
	out_node.angst_name = NODE_N[node_name].angst_name

	#List of nodes included in this clade
	read = []
	dec = NODE_N[node_name].dec
	dec_angst = {}
	dec_angst[NODE_N[node_name].name] = ""
	for j in dec:
		if ancestral_flag == 'an':
			if NODE_N[j].genome == 0:
				read.append(j)
				dec_angst[NODE_N[j].name] = ""
				out_node.include[j] = ''
			read.append(node_name)
		elif ancestral_flag == 'g':
			if NODE_N[j].genome == 1:
				read.append(j)
				dec_angst[NODE_N[j].name] = ""
				out_node.include[j] = ''
		elif ancestral_flag == 'all':
			read.append(j)
			dec_angst[NODE_N[j].name] = ""
			out_node.include[j] = ''
			read.append(node_name)
		else:
			print 'Type flag must be an (ancestral nodes), g (genome nodes), or all (all nodes)'

	#For each node in the subclade
	for j in read:
		out_node.los_total += NODE_N[j].los_total
		out_node.dup_total += NODE_N[j].dup_total
		out_node.brn_total += NODE_N[j].brn_total
		out_node.spc_total += NODE_N[j].spc_total
		out_node.event_total += NODE_N[j].event_total
		out_node.total += NODE_N[j].total

		out_node.lgt_total_source += NODE_N[j].lgt_total_source
		out_node.lgt_total_target += NODE_N[j].lgt_total_target

		for cat in NODE_N[j].cats:
			if not out_node.cats.has_key(cat):
				out_node.cats[cat] = 0
			out_node.cats[cat] += NODE_N[j].cats[cat]

		for g in NODE_N[j].lgt_source:
			if not out_node.lgt_source.has_key(g):
				out_node.lgt_source[g] = []
			out_node.lgt_source[g] += NODE_N[j].lgt_source[g]

		for g in NODE_N[j].lgt_target:
			if not out_node.lgt_target.has_key(g):
				out_node.lgt_target[g] = []
			out_node.lgt_target[g] += NODE_N[j].lgt_target[g]

		for g in NODE_N[j].lgt_source_gene:
			if not out_node.lgt_source_gene.has_key(g):
				out_node.lgt_source_gene[g] = 0
			out_node.lgt_source_gene[g] += NODE_N[j].lgt_source_gene[g]

		for g in NODE_N[j].lgt_target_gene:
			if not out_node.lgt_target_gene.has_key(g):
				out_node.lgt_target_gene[g] = 0
			out_node.lgt_target_gene[g] += NODE_N[j].lgt_target_gene[g]
	return out_node


#
#Rates of LGT between nodes (modify as needed)
#
def calc_ratesA():
	"""Calculating rates of LGT between clades, needs to be modified for correct boundaries"""
	global NODE_N

	AA_perc = []	#LGT percent from A
	AB_perc = []	#LGT percent from B
	AA = []
	AB = []
	ABasal_perc = []	#LGT percent from basal strepts
	ABasal = []
	AActino = []
	AActino_perc = []
	A_strept = []
	total_lgt_A = []
	gCount = 0
	AGE = {}
	AGE['A'] = 120.273			#Multilocus: 120.273	Astral: 125.4957
	AGE['B'] = 118.25			#Multilocus: 118.25		Astral: 126.5378
	AGE['AB'] = 162.1488		#Multilocus: 162.1488	Astral: 141.3703
	AGE['basal'] = 333.4095		#Multilocus: 333.4095	Astral: 321.6617
	AGE['actino'] = 328.7652		#Multilcous: 328.7652	Astral: 369.024

	b_node = 'NODE89'			#Multilocus: 89		Astral: 126
	a_node = 'NODE126'			#Multilocus: 126	Astral: 77
	strept_node = 'NODE133'		#Multilocus: 133	Astral: 132

	for n in NODE_N:
		total = 0
		A = 0
		B = 0
		basal = 0
		strept = 0
		actino = 0
		if NODE_N[a_node].dec.has_key(n):
			if NODE_N[n].genome == 1:
				gCount += 1

			for lgt in NODE_N[n].lgt_target:
				if NODE_N[b_node].dec.has_key(lgt):
					B += float(len(NODE_N[n].lgt_target[lgt]))
					strept += float(len(NODE_N[n].lgt_target[lgt]))
				elif NODE_N[a_node].dec.has_key(lgt):
					A += float(len(NODE_N[n].lgt_target[lgt]))
					strept += float(len(NODE_N[n].lgt_target[lgt]))
				elif NODE_N[strept_node].dec.has_key(lgt):
					basal += float(len(NODE_N[n].lgt_target[lgt]))
					strept += float(len(NODE_N[n].lgt_target[lgt]))
				else:
					actino += float(len(NODE_N[n].lgt_target[lgt]))
				total += float(len(NODE_N[n].lgt_target[lgt]))
		if total > 0:
			total_lgt_A.append(total)
			AA_perc.append(A / total)
			AA.append(A)
			AB_perc.append(B / total)
			AB.append(B)
			ABasal_perc.append(basal / total)
			ABasal.append(basal)
			A_strept.append(strept / total)
			AActino_perc.append(actino / total)
			AActino.append(actino)

	print 'Total from A: {0}'.format(sum(AA))
	print 'Total from B: {0}'.format(sum(AB))
	print 'Total from Basal: {0}'.format(sum(ABasal))
	print 'Total from Actino: {0}'.format(sum(AActino))
	print 'Total % from A: {0}'.format(sum(AA) / sum(total_lgt_A))
	print 'Total % from B: {0}'.format(sum(AB) / sum(total_lgt_A))
	print 'Total % from Basal: {0}'.format(sum(ABasal) / sum(total_lgt_A))
	print 'Total % from Actino: {0}\n'.format(sum(AActino) / sum(total_lgt_A))

	perNodeAA = sum(AA) / len(AA)
	perNode_ageAA = perNodeAA / AGE['A']
	print 'Per-node from A: {0}, {1} per my'.format(perNodeAA,perNode_ageAA)

	perNodeAB = sum(AB) / len(AA)
	perNode_ageAB = perNodeAB / AGE['AB']
	print 'Per-node from B: {0}, {1} per my'.format(perNodeAB,perNode_ageAB)

	perNodeABasal = sum(ABasal) / len(AA)
	perNode_ageABasal = perNodeABasal / AGE['basal']
	print 'Per-node from Basal: {0}, {1} per my'.format(perNodeABasal,perNode_ageABasal)

	perNodeAActino = sum(AActino) / len(AA)
	perNode_ageAActino = perNodeAActino / AGE['actino']
	print 'Per-node from Actino: {0}, {1} per my\n'.format(perNodeAActino,perNode_ageAActino)


	perGenomeAA = sum(AA) / gCount
	perGenome_ageAA = perGenomeAA / AGE['A']
	print 'Per-genome from A: {0}, {1} per my'.format(perGenomeAA,perGenome_ageAA)

	perGenomeAB = sum(AB) / gCount
	perGenome_ageAB = perGenomeAB / AGE['AB']
	print 'Per-genome from B: {0}, {1} per my'.format(perGenomeAB,perGenome_ageAB)

	perGenomeABasal = sum(ABasal) / gCount
	perGenome_ageABasal = perGenomeABasal / AGE['basal']
	print 'Per-genome from Basal: {0}, {1} per my'.format(perGenomeABasal,perGenome_ageABasal)

	perGenomeAActino = sum(AActino) / gCount
	perGenome_ageAActino = perGenomeAActino / AGE['actino']
	print 'Per-genome from Actino: {0}, {1} per my'.format(perGenomeAActino,perGenome_ageAActino)


def calc_ratesB():
	"""Calculating rates of LGT between clades, needs to be modified for correct boundaries"""
	global NODE_N

	BA_perc = []	#LGT percent from A
	BB_perc = []	#LGT percent from B
	BA = []
	BB = []
	BBasal_perc = []	#LGT percent from basal strepts
	BBasal = []
	BActino = []
	BActino_perc = []
	B_strept = []
	total_lgt_B = []
	gCount = 0
	AGE = {}
	AGE['A'] = 120.273			#Multilocus: 120.273	Astral: 125.4957
	AGE['B'] = 118.25			#Multilocus: 118.25		Astral: 126.5378
	AGE['AB'] = 162.1488		#Multilocus: 162.1488	Astral: 141.3703
	AGE['basal'] = 333.4095		#Multilocus: 333.4095	Astral: 321.6617
	AGE['actino'] = 328.7652		#Multilcous: 328.7652	Astral: 369.024

	b_node = 'NODE89'			#Multilocus: 89		Astral: 126
	a_node = 'NODE126'			#Multilocus: 126	Astral: 77
	strept_node = 'NODE133'		#Multilocus: 133	Astral: 132

	for n in NODE_N:
		total = 0
		A = 0
		B = 0
		basal = 0
		strept = 0
		actino = 0
		if NODE_N[b_node].dec.has_key(n):
			if NODE_N[n].genome == 1:
				gCount += 1

			for lgt in NODE_N[n].lgt_target:
				if NODE_N[b_node].dec.has_key(lgt):
					B += float(len(NODE_N[n].lgt_target[lgt]))
					strept += float(len(NODE_N[n].lgt_target[lgt]))
				elif NODE_N[a_node].dec.has_key(lgt):
					A += float(len(NODE_N[n].lgt_target[lgt]))
					strept += float(len(NODE_N[n].lgt_target[lgt]))
				elif NODE_N[strept_node].dec.has_key(lgt):
					basal += float(len(NODE_N[n].lgt_target[lgt]))
					strept += float(len(NODE_N[n].lgt_target[lgt]))
				else:
					actino += float(len(NODE_N[n].lgt_target[lgt]))
				total += float(len(NODE_N[n].lgt_target[lgt]))

		if total > 0:
			total_lgt_B.append(total)
			BA_perc.append(A / total)
			BA.append(A)
			BB_perc.append(B / total)
			BB.append(B)
			BBasal_perc.append(basal / total)
			BBasal.append(basal)
			B_strept.append(strept / total)
			BActino_perc.append(actino / total)
			BActino.append(actino)

	print 'Total from A: {0}'.format(sum(BA))
	print 'Total from B: {0}'.format(sum(BB))
	print 'Total from Basal: {0}'.format(sum(BBasal))
	print 'Total from Actino: {0}'.format(sum(BActino))
	print 'Total % from A: {0}'.format(sum(BA) / sum(total_lgt_B))
	print 'Total % from B: {0}'.format(sum(BB) / sum(total_lgt_B))
	print 'Total % from Basal: {0}'.format(sum(BBasal) / sum(total_lgt_B))
	print 'Total % from Actino: {0}\n'.format(sum(BActino) / sum(total_lgt_B))

	perNodeBA = sum(BA) / len(BA)
	perNode_ageBA = perNodeBA / AGE['AB']
	print 'Per-node from A: {0}, {1} per my'.format(perNodeBA,perNode_ageBA)

	perNodeBB = sum(BB) / len(BB)
	perNode_ageBB = perNodeBB / AGE['B']
	print 'Per-node from B: {0}, {1} per my'.format(perNodeBB,perNode_ageBB)

	perNodeBBasal = sum(BBasal) / len(BA)
	perNode_ageBBasal = perNodeBBasal / AGE['basal']
	print 'Per-node from Basal: {0}, {1} per my'.format(perNodeBBasal,perNode_ageBBasal)

	perNodeBActino = sum(BActino) / len(BA)
	perNode_ageBActino = perNodeBActino / AGE['actino']
	print 'Per-node from Actino: {0}, {1} per my\n'.format(perNodeBActino,perNode_ageBActino)


	perGenomeBA = sum(BA) / gCount
	perGenome_ageBA = perGenomeBA / AGE['AB']
	print 'Per-genome from A: {0}, {1} per my'.format(perGenomeBA,perGenome_ageBA)

	perGenomeBB = sum(BB) / gCount
	perGenome_ageBB = perGenomeBB / AGE['B']
	print 'Per-genome from B: {0}, {1} per my'.format(perGenomeBB,perGenome_ageBB)

	perGenomeBBasal = sum(BBasal) / gCount
	perGenome_ageBBasal = perGenomeBBasal / AGE['basal']
	print 'Per-genome from Basal: {0}, {1} per my'.format(perGenomeBBasal,perGenome_ageBBasal)

	perGenomeBActino = sum(BActino) / gCount
	perGenome_ageBActino = perGenomeBActino / AGE['actino']
	print 'Per-genome from Actino: {0}, {1} per my'.format(perGenomeBActino,perGenome_ageBActino)

def geneLGTDist(minDist,conv,outputFile,print_fams='ALL'):
	"""Identifies recent LGT events: genes that were transfered to their current genome 
or a recent ancestor within a certain evolutionary distance."""
	global LEAF
	global NODE_N
	global FAMS
	global tree

	REC = {}			#List of genes with recent LGTs
	LGTDist = {}		#Dict of genes in each genome, whether they have recent LGTs or not

	class LGT_Genome():
		def __init__(self):
			self.name = ''
			self.lgt = {}
			self.dup = {}
			self.spc = {}
			self.total = {}
			self.los = {}
			self.old_lgt = {}
			self.extant_lgt = {}
			self.old_los = {}
			self.extant_los = {}
			self.genes = {} 

			for f in FAMS:
				self.lgt[f] = []
				self.dup[f] = []
				self.spc[f] = []
				self.total[f] = []
				self.los[f] = 0
				self.old_lgt[f] = []
				self.extant_lgt[f] = []
				self.old_los[f] = []
				self.extant_los[f] = []

				# Records most recent event:
				# LGT_GENOME[genome].genes[fam][gene] = [event_type,age]
				self.genes[f] = {}

	genomes = []
	for n in NODE_N:
		if NODE_N[n].genome == 1:
			REC[n] = LGT_Genome()
			REC[n].name = n
			genomes.append(NODE_N[n].name)

	if print_fams:
		RUN_FAMS = {}
		if print_fams == 'ALL':
			for fam in FAMS:
				RUN_FAMS[fam] = ''
		else:
			for i in print_fams:
				if FAMS.has_key(i):
					RUN_FAMS[i] = ''

	for genome in LEAF:
		LGTDist[genome] = {}
		for g in LEAF[genome]:
			fam = LEAF[genome][g].fam

			#In case some fams werent run through angst
			if not RUN_FAMS.has_key(fam):
				continue

			REC[genome].total[fam].append(g)
			REC[genome].genes[fam][g] = ['spc',minDist]

			event = 0
			for t in LEAF[genome][g].lgt:
				target_dist = tree.distance(genome,t[1]) * conv		#Age from LGT event to extant gene
				if target_dist <= minDist:
					LGTDist[genome][g] = target_dist
					REC[genome].lgt[fam].append(LEAF[genome][g].name)
					event = 1
					if t[1] != genome:
						REC[genome].old_lgt[fam].append(LEAF[genome][g].name)
					else:
						REC[genome].extant_lgt[fam].append(LEAF[genome][g].name)

					if target_dist <= REC[genome].genes[fam][g][1]:
						REC[genome].genes[fam][g] = ['lgt',target_dist,t[0]]


			for t in LEAF[genome][g].dup:
				target_dist = tree.distance(genome,t) * conv		#Age from LGT event to extant gene
				if target_dist <= minDist:
					REC[genome].dup[fam].append(LEAF[genome][g].name)

					if target_dist <= REC[genome].genes[fam][g][1]:
						REC[genome].genes[fam][g] = ['dup',target_dist]
					event = 1

			if event == 0:
				REC[genome].spc[fam].append(LEAF[genome][g].name)

	for genome in genomes:
		path = tree.trace(tree.root,genome)
		for node in path:
			if (tree.distance(genome,node.name) * conv) <= minDist:
				for fam in RUN_FAMS:
					REC[genome].los[fam] += NODE_N[node.name].los[fam]

	os.system('mkdir -p {0}'.format(outputFile))
	out_stats = open('{1}/fam_stats_{0}.tbl'.format(minDist,outputFile),'w')
	out_stats.write('Fam\tLGT\tDup\tSpc\tLos\tTotal\n')
	out_ratios = open('{1}/fam_ratios_{0}.tbl'.format(minDist,outputFile),'w')
	out_ratios.write('Fam\tLGT_per_gene\tDup_per_gene\tSpc_per_gene\tLos_per_gene\tTotal\n')
	if print_fams:

		for fam in RUN_FAMS:
			lgt_total = 0.0
			dup_total = 0.0
			los_total = 0.0
			spc_total = 0.0
			count_total = 0.0
			out = open('{2}/events_{1}my_{0}.tbl'.format(fam,minDist,outputFile),'w')
			out.write('Genome\tLGT\tDup\tSpc\tLos\tTotal\n')
			for g in REC:
				lgt_g = 0
				dup_g = 0
				los_g = REC[g].los[fam]
				spc_g = 0
				count_g = len(REC[g].genes[fam])
				count_total += len(REC[g].genes[fam])
				los_total += REC[g].los[fam]

				for gene in REC[g].genes[fam]:
					if REC[g].genes[fam][gene][0] == 'lgt':
						lgt_total += 1
						lgt_g += 1
					elif REC[g].genes[fam][gene][0] == 'dup':
						dup_total += 1
						dup_g += 1
					elif REC[g].genes[fam][gene][0] == 'spc':
						spc_total += 1
						spc_g += 1

				data = [g,lgt_g,dup_g,spc_g,los_g,count_g]
				data2 = [str(x) for x in data]
				out.write('\t'.join(data2) + '\n')
			out.close()

			data = [lgt_total,dup_total,spc_total,los_total,count_total]
			data2 = [str(x) for x in data]

			ratios = [x / count_total for x in data]
			ratios2 = [str(x) for x in ratios]

			out_stats.write(fam + '\t' + '\t'.join(data2) + '\n')
			out_ratios.write(fam + '\t' + '\t'.join(ratios2) + '\n')
	out_stats.close()

	return (REC,LGTDist)




def clusterLGT(inFolder,minDist,conv,outFolder):
	"""Identifies percent of genes in each gene cluster that have been acquired within a certain timeframe"""
	#inFolder should be a folder of gene lists
	import collections as cl

	global LEAF
	global NODE_N
	global FAMS
	global tree

	class Cluster():
		def __init__(self):
			self.name = ''
			self.genes = {}
			self.lgt = {}
			self.dup = {}

	GENES = {}
	#Need to swap keys for values in the FAMS.genes dictionary for gene lookup
	for fam in FAMS:
		for gene in FAMS[fam].genes:
			GENES[gene] = fam

	#Getting the genes in each cluster
	print '\nReading genes in clusters...'
	CLUST = {}
	clusterFiles = os.listdir(inFolder)
	for f in clusterFiles:
		name = f.split('.')[0]
		CLUST[name] = Cluster()
		CLUST[name].name = name

		with open('{0}/{1}'.format(inFolder,f),'r') as c_data:
			for ln in c_data:
				ln = ln.translate(None,'>\n')
				ln = ln.split('|')[-1]
				CLUST[name].genes[ln] = ''

	#Only look at LGT events for the fams covered by genes in the clusters,
	#speeds up analysis for AnGST runs with high numbers of gene families
	search_fams = cl.deque()
	for c in CLUST:
		for g in CLUST[c].genes:
			if GENES.has_key(g):
				search_fams.append(GENES[g])
	search_fams = set(search_fams)

	#Calculating LGT events using the function above
	print 'Calculating LGT events...'
	rawGeneOut = 'raw_' + outFolder
	REC, LGTDist = geneLGTDist(minDist,conv,rawGeneOut,print_fams=search_fams)

	print'Identifying recent LGTs affecting cluster genes...'
	GENOME = {}
	for c in CLUST:
		for gene in CLUST[c].genes:
			if GENES.has_key(gene):
				fam = GENES[gene]
				genome = gene.split('_')[0]
				if not GENOME.has_key(genome):
					GENOME[genome] = {}
				GENOME[genome][c] = CLUST[c]

				#Looking for each gene in the cluster in the 'recently lgted' list
				#Messy way to do this, but much faster than 'is x in list'
				GENE_TRACKER = {}
				for thing in REC[genome].lgt[fam]:
						if REC[genome].genes[fam][thing][0] == 'lgt':
							GENE_TRACKER[thing] = REC[genome].genes[fam][thing][2]
				if GENE_TRACKER.has_key(gene):
					GENOME[genome][c].lgt[gene] = GENE_TRACKER[gene]

	os.system('mkdir -p {0}'.format(outFolder))
	out = open('{0}/00_column_labels.txt'.format(outFolder),'w')
	header = ('Cluster','Total_LGT_Genes','Ratio_lgt','Primary_source','LGT_genes_from_primary',
		'Ratio_primary_lgt','Total_Genes')
	out.write('\t'.join(header))
	os.system('mkdir -p {0}/'.format(outFolder))
	for g in GENOME:
		out = open('{0}/{1}.cluster_lgt'.format(outFolder,g),'w')
		for c in GENOME[g]:
			lgtCount = 0.0
			totalGene = 0.0
			SOURCE = {}
			for gene in GENOME[g][c].genes:
				if GENES.has_key(gene):
					genome = gene.split('_')[0]
					if GENOME[g][c].lgt.has_key(gene):
						lgtCount += 1
						lgtSource = GENOME[g][c].lgt[gene]
						if not SOURCE.has_key(lgtSource):
							SOURCE[lgtSource] = 0.0
						SOURCE[lgtSource] += 1
					totalGene += 1
			if totalGene:
				lgtRatio = lgtCount / totalGene
				mainSource = ''
				maxLgt = 0
				for i in SOURCE:
					if SOURCE[i] > maxLgt:
						mainSource = i
						maxLgt = SOURCE[i]
				primaryRatio = maxLgt / totalGene
				data = (c,lgtCount,lgtRatio,mainSource,maxLgt,primaryRatio,totalGene)
				data = [str(x) for x in data]

				out.write('\t'.join(data) + '\n')
	out.close()

	return CLUST,REC,LGTDist




def tip_lgt_rate(NODE_N,clade=''):
	lgt_tip_perc = []
	lgt_in_perc = []
	lgt_tip = []
	lgt_in = []

	for n in NODE_N:
		if clade and not NODE_N[clade].dec.has_key(n):
			continue

		total = NODE_N[n].total
		lgt = 0.0
		for i in NODE_N[n].lgt_target:
			lgt += len(NODE_N[n].lgt_target[i])

		if 'NODE' not in n:
			lgt_tip_perc.append(lgt/total)
			lgt_tip.append(lgt)
		else:
			lgt_in_perc.append(lgt/total)
			lgt_in.append(lgt)

		var_perc = stats.mannwhitneyu(lgt_tip_perc,lgt_in_perc)
		mean_in_perc = np.mean(lgt_in_perc)
		mean_tip_perc = np.mean(lgt_tip_perc)
		std_tip_perc = np.std(lgt_tip_perc)
		std_in_perc = np.std(lgt_in_perc)

		var_num = stats.mannwhitneyu(lgt_tip,lgt_in)
		mean_in = np.mean(lgt_in)
		mean_tip = np.mean(lgt_tip)
		std_tip = np.std(lgt_tip)
		std_in = np.std(lgt_in)

	print '\nTip percent LGT:\t{0} +/- {1}'.format(mean_tip_perc,std_tip_perc)
	print 'Internal percent LGT:\t{0} +/- {1}'.format(mean_in_perc,std_in_perc)
	print 'Pval percent:\t\t{0}\n'.format(var_perc[1])
	print 'Tip genes LGT:\t\t{0} +/- {1}'.format(mean_tip,std_tip)
	print 'Internal genes LGT:\t{0} +/- {1}'.format(mean_in,std_in)
	print 'Pval genes:\t\t{0}'.format(var_num[1])

def fet_clade_lgt(node):
	"""FET for correlation between LGTs into a clade and the source of those LGTs"""
	global NODE_N
	
	#		LGT_Recipient_Clade
	#					node	!node
	#LGT		node  	a		  b
	#Source		!node  	c		  d
	#clade

	DEC = NODE_N[node].dec

	a = 0
	b = 0
	c = 0
	d = 0
	for n in NODE_N:
		for i in NODE_N[n].lgt_target:

			#A few excess 'in' or 'not in' here, but more explicit
			if i in DEC and n in DEC:
				a += 1
			elif i in DEC and n not in DEC:
				c += 1
			elif i not in DEC and n in DEC:
				b += 1
			elif i not in DEC and n not in DEC:
				d += 1
	odds, pval = stats.fisher_exact([[a, b], [c, d]])
	return odds,pval,[a,b,c,d]

def fet_leaf_clade_lgt(node,geneList):
	"""FET for correlation between LGT events into a clade vs the source,
using only the LGT events in the histroy of a specific set of genes"""
	global NODE_N
	global LEAF
	
	#		LGT_Recipient_Clade
	#					node	!node
	#LGT		node  	a		  b
	#Source		!node  	c		  d
	#clade

	genes = []
	with open(geneList,'r') as f:
		for ln in f:
			ln = ln.translate(None,'\n>')
			genes.append(ln.split('|')[-1])

	DEC = NODE_N[node].dec
	TERMINAL = NODE_N[node].term_list

	a = 0
	b = 0
	c = 0
	d = 0
	for g in genes:
		genome = g.split('_')[0]
		if LEAF[genome].has_key(g):
			for i in LEAF[genome][g].lgt:
				source = i[0]
				target = i[1]

				#A few excess 'in' or 'not in' here, but more explicit
				if source in DEC and target in DEC:
					a += 1
				elif source in DEC and target not in DEC:
					c += 1
				elif source not in DEC and target in DEC:
					b += 1
				elif source not in DEC and target not in DEC:
					d += 1
	odds, pval = stats.fisher_exact([[a, b], [c, d]])
	return odds,pval,[a,b,c,d]



def tree_networks(NODE_N,angst_folder):
	"""Generates a composite network for each gene family, the species tree with additional zero-length edges
	for LGT events. Enables tracing of the history of extant genes through the phylogenetic network  using nx.dijkstra_path_length"""
	import networkx as nx

	coreTree = Phylo.read('{0}/species_tree_labels.newick'.format(angst_folder),'newick')
	coreNet = nx.Graph()
	for i in coreTree.get_terminals():
		coreNet.add_node(i.name)
	for i in coreTree.get_nonterminals():
		coreNet.add_node(i.name)
	for nodeA in coreTree.get_nonterminals():
		for nodeB in coreTree.get_nonterminals():
			if nodeA.is_parent_of(nodeB) and len(coreTree.trace(nodeA,nodeB)) == 2:
				dist = coreTree.distance(nodeA,nodeB)
				coreNet.add_edge(nodeA.name,nodeB.name,weight=dist)
		for nodeB in coreTree.get_terminals():
			if nodeA.is_parent_of(nodeB):
				dist = coreTree.distance(nodeA,nodeB)
				coreNet.add_edge(nodeA.name,nodeB.name,weight=dist)