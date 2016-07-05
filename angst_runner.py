import sys
import os
import re
import numpy as np
from multiprocessing import Pool
from ete2 import Tree
from Bio import Phylo

usage = """\nRuns angst on a set of gene trees\n
<species tree>
<folder of gene trees>
<output folder>
<cores>
<format trees? (Y/N)>
<ultrametric tree? (Y/N)>
<LGT Penalty (default 3)>
<Duplication Penalty (default 2)>
<Timeout max ram (in gb), or 0 if not using timeout>\n\n"""
argv = sys.argv[1:]

if len(argv) == 9:
	species_tree = argv.pop(0)
	gene_tree_folder = argv.pop(0)
	output_folder = argv.pop(0)
	cores = int(argv.pop(0))
	format_spec_flag = argv.pop(0)
	ultra_flag = argv.pop(0)
	lgt_pen = float(argv.pop(0))
	dup_pen = float(argv.pop(0))
	timeout = int(float(argv.pop(0)) * 1048576)		#Conversion to kb, which is what timeout needs
	max_ram_mb = float(timeout) / 1000	#Conversion to MB for ram cutoff
else:
	sys.exit(usage)

#
#Install location for AnGST primary script
#

if timeout > 0:
	timeout_install = "~/tools_and_software/timeout/timeout -m {0}".format(timeout)

angst_install = "~/tools_and_software/angst_slots_cython/angst_lib/AnGST.py"

#Formula to approximate ram usage for a gene family.
def approx_ram(gene_count):
	return np.exp((1.944377 * np.log(float(gene_count))) - 2.435937)


#Runs commands for building directory and running AnGST for each gene tree
def runner(gene):
	name = gene
	name = re.sub("^\S+\.", "", name)

	ram_usage = approx_ram(COUNT[gene]) * 0.8
	if ram_usage > max_ram_mb and timeout > 0:
		print 'Skipping {0} - predicted ram usage: {1} genes, {2}mb ram'.format(gene,COUNT[gene],ram_usage)
		return

	#Writing angst input file
	oIN = open("{0}/{1}_master_input".format(output_folder,name), 'w')
	oIN.write("species={0}\n".format(species_tree))
	oIN.write("gene={0}/{1}\n".format(gene_tree_folder,gene))
	oIN.write("output={0}/{1}_angst/\n".format(output_folder,name))
	oIN.write("penalties={0}/penalties".format(output_folder,name))
	if ultra_flag == "Y":
		oIN.write("\nultrametric=True".format(output_folder,name))
	oIN.close()

	if timeout > 0:
		os.system("{3} python {0} {1}/{2}_master_input".format(angst_install, output_folder,name,timeout_install))
	else:
		os.system("python {0} {1}/{2}_master_input".format(angst_install, output_folder,name))
	os.system("mv {0}/{1}_master_input {0}/master_inputs/{1}_master_input".format(output_folder,name))

#
#Script Start
#
os.system("mkdir -p {0}".format(output_folder))
os.system('mkdir -p {0}/master_inputs/'.format(output_folder))

#Writing penalties file
PEN = open("{0}/penalties".format(output_folder), 'w')
PEN.write("hgt: {0}\ndup: {1}\nlos: 1.0\nspc: 0.0\n".format(lgt_pen,dup_pen))
PEN.close()

#Reformatting species tree if required
if "Y" in format_spec_flag or "y" in format_spec_flag:
	bootstrap = re.compile("\[\d+\]")
	sp_tree = open(species_tree, 'r')
	sp_tree = sp_tree.readlines()
	sp_tree = ''.join(sp_tree)
	sp_tree = sp_tree.translate(None, "\n\r")
	sp_tree = sp_tree.replace(";", "")

	sp_tree = sp_tree.replace("_", "/")
	sp_tree = bootstrap.sub("",sp_tree)
	sp_clean = open("{0}/species_tree_clean".format(output_folder), 'w')
	sp_clean.write("({0}:0.0);".format(sp_tree))
	sp_clean.close()
	species_tree = "{0}/species_tree_clean".format(output_folder)
elif "N" in format_spec_flag or "n" in format_spec_flag:
	os.system("cp {0} {1}/species_tree_clean".format(species_tree, output_folder))
	species_tree = "{0}/species_tree_clean".format(output_folder)
else:
	species_tree = "{0}/species_tree_clean".format(output_folder)

#
#Reformatting gene trees and counting the number of taxa in each tree
#
os.system("mkdir -p {0}/gene_trees_rough".format(output_folder))
os.system('mkdir -p {0}/gene_trees'.format(output_folder))
gene_list = os.listdir(gene_tree_folder)
COUNT = {}

if "Y" in format_spec_flag or "y" in format_spec_flag:
	#Bio.Phylo: rename genes, removing |s and plasmid names
	for i in gene_list:
		tList = []
		for tree in Phylo.parse("{0}/{1}".format(gene_tree_folder,i),'newick'):
			tList.append(tree)
		for t in tList:

			#formatting names
			for node in t.get_terminals():
				name = node.name.split('|')
				if len(name) == 2:
					node.name = node.name.split('|')[1]
				elif len(name) == 3:
					node.name = node.name.split('|')[2]
		Phylo.write(tList,"{0}/gene_trees_rough/{1}".format(output_folder, str(i)),'newick')

	#ete2: remove polytomies, unroot, and format for AnGST's parser to read
	gene_trees = os.listdir('{0}/gene_trees_rough'.format(output_folder))
	for gene_tree in gene_trees:
		with open('{0}/gene_trees_rough/{1}'.format(output_folder,gene_tree),'r') as f:
			out = open('{0}/gene_trees/{1}'.format(output_folder,gene_tree),'w')
			max_genes = 0

			for ln in f:
				polyFlag = 0
				ln = ln.replace("\n",'')
				t = Tree(ln)
				if len(t.get_leaves()) > max_genes:
					max_genes = len(t.get_leaves())
				if max_genes > 2:
					t.resolve_polytomy(recursive=True)
					t.unroot()
					out.write('{0}\n'.format(t.write(format=5)))
			out.close()
		if max_genes < 3:
			os.system('rm {0}/gene_trees/{1}'.format(output_folder,gene_tree))


gene_tree_folder = '{0}/gene_trees'.format(output_folder)
gene_list = os.listdir(gene_tree_folder)
for i in gene_list:
	for tree in Phylo.parse("{0}/{1}".format(gene_tree_folder,i),'newick'):
		if not COUNT.has_key(i):
			COUNT[i] = len(tree.get_terminals())
			break

#Running AnGST
gene_list = os.listdir(gene_tree_folder)
p=Pool(cores)
for gene in gene_list:
	p.apply_async(runner, args=(gene,))
	#pass
p.close()
p.join()

#
#Compiling data from AnGST run
#
angst_data = os.listdir(output_folder)
os.system("mkdir -p {0}/data_files".format(output_folder))
out_sum = open("{0}/00_data_summary.out".format(output_folder),'w')
summary = ["hgt inferred:", "dup inferred:", "los inferred:", "spc inferred:", "gene count:",
"run-time", "last reconciliation"]

for fi in angst_data:
	if "_angst" in fi:
		gene_data = os.listdir("{0}/{1}".format(output_folder,fi))
		gene = fi
		gene = gene.replace("_angst", "")
		out = open("{0}/data_files/{1}_data".format(output_folder,gene), 'w')
		out_sum.write("####    {0}    ####\n".format(gene))

		for f in gene_data:
			stats_flag = 0
			name = str(f)
			name = name.replace("AnGST.", "")
			if "stats" in name:
				stats_flag = 1

			out.write("\n#########\t{0}\t#########\n".format(name))
			data = open("{0}/{1}/{2}".format(output_folder,fi,f),'r')

			data = data.readlines()
			for i in data:
				if stats_flag != 1:
					i = i.replace("/","_")
				out.write(i)
				if stats_flag == 1:
					for j in summary:
						if j in i:
							out_sum.write(i)
			out.write("\n############################\n\n")
		out.close()
out_sum.close()

os.system("rm -f -r {0}/*_angst".format(output_folder))
os.system("rm -f -r {0}/penalties".format(output_folder))

