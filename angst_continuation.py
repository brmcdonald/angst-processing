import sys,os,re
import numpy as np
from Bio import Phylo
from multiprocessing import Process, Queue, Pool

usage = """\nReads an existing AnGST run and finishes any families not already run\n
<Angst folder>
<cores>
<ultrametric tree? (Y/N)>
<LGT Penalty (default 3)>
<Duplication Penalty (default 2)>
<Timeout max ram (in gb), or 0 if not using timeout>
<Optional: priority fams, separated by commas>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(argv) == 6:
	angst_folder = argv.pop(0)
	cores = int(argv.pop(0))
	ultra_flag = argv.pop(0)
	lgt_pen = float(argv.pop(0))
	dup_pen = float(argv.pop(0))
	timeout = int(float(argv.pop(0)) * 1048576) #Conversion to kb, which is what timeout needs
	max_ram_mb = float(timeout) / 1000			#Conversion to MB for ram cutoff
	runFirst = []
elif len(argv) == 7:
	angst_folder = argv.pop(0)
	cores = int(argv.pop(0))
	ultra_flag = argv.pop(0)
	lgt_pen = float(argv.pop(0))
	dup_pen = float(argv.pop(0))
	timeout = int(float(argv.pop(0)) * 1048576) #Conversion to kb, which is what timeout needs
	max_ram_mb = float(timeout) / 1000	#Conversion to MB for ram cutoff
	runFirst = argv.pop(0).split(',')
else:
	sys.exit(usage)

#
#Install location for AnGST primary script
#
if timeout > 0:
	timeout_install = "~/tools_and_software/timeout/timeout -m {0}".format(timeout)

angst_install = "~/tools_and_software/angst_slots_cython/angst_lib/AnGST.py"

def approx_ram(gene_count):
	return np.exp((1.944377 * np.log(float(gene_count))) - 2.435937)

#Runs commands for building directory and running AnGST for each gene tree
def runner(gene):
	name = gene
	name = re.sub("^\S+\.", "", name)

	ram_usage = approx_ram(COUNT[gene]) * 0.8
	if ram_usage > max_ram_mb and timeout > 0:
		print '\nSkipping {0} - predicted ram usage: {1} genes, {2}mb ram'.format(gene,COUNT[gene],ram_usage)
		return
	else:
		print '\nRunning {0}'.format(gene)

	#Writing angst input file
	oIN = open("{0}/{1}_master_input".format(angst_folder,name), 'w')
	oIN.write("species={0}\n".format(species_tree))
	oIN.write("gene={0}/gene_trees/{1}\n".format(angst_folder,gene))
	oIN.write("output={0}/{1}_angst/\n".format(angst_folder,name))
	oIN.write("penalties={0}/penalties".format(angst_folder,name))
	if ultra_flag == "Y":
		oIN.write("\nultrametric=True".format(angst_folder,name))
	oIN.close()

	if timeout > 0:
		os.system("{3} python {0} {1}/{2}_master_input".format(angst_install, angst_folder,name,timeout_install))
	else:
		os.system("python {0} {1}/{2}_master_input".format(angst_install, angst_folder,name))
	os.system("mv {0}/{1}_master_input {0}/master_inputs/{1}_master_input".format(output_folder,name))

#
#Script Start
#

#Writing penalties file
PEN = open("{0}/penalties".format(angst_folder), 'w')
PEN.write("hgt: {0}\ndup: {1}\nlos: 1.0\nspc: 0.0\n".format(lgt_pen,dup_pen))
PEN.close()
species_tree = '{0}/species_tree_clean'.format(angst_folder)

#
#Reading gene tree list
#
GENES = {}
gene_trees = os.listdir('{0}/gene_trees'.format(angst_folder))
for f in gene_trees:
	f_clean = re.sub('_boot_\d+\s*','',f)
	GENES[f_clean] = f
#
#Reading finished gene families
#
GENES_DONE = {}
folder_data = os.listdir(angst_folder)

#If the run completely finished, read the data summary
if '00_data_summary.out' in folder_data:
	data = open('{0}/00_data_summary.out'.format(angst_folder),'r')
	data = data.readlines()
	for ln in data:
		if'###' in ln:
			m = re.search('\#\#\#\s+(\S+)\s+\#\#\#',ln)
			if m:
				gene = m.group(1)
				g_clean = re.sub('_boot_\d+','',gene)
				GENES_DONE[g_clean] = ''
			else:
				print 'Cant find gene fam in {0}'.format(ln.replace('\n',''))

#If things aren't in the data summary or it didnt finish all the way, read the output
#folders
for i in folder_data:
	if '_angst' in i and '_boot_' in i:
		i = re.sub('_boot_\d+_angst','',i)
		GENES_DONE[i] = ''

#
#Figuring out which gene fams are missing
#

COUNT = {}
gene_list = []
priority_gene_list = []
for i in GENES:
	if not GENES_DONE.has_key(i):
		gene_list.append(GENES[i])

		for tree in Phylo.parse('{0}/gene_trees/{1}'.format(angst_folder,GENES[i]),'newick'):
			if not COUNT.has_key(GENES[i]):
				COUNT[GENES[i]] = len(tree.get_terminals())
				break

print '\n***\t{0} gene families not run yet\t***\n'.format(len(gene_list))

if runFirst:
	remove = []
	for gene in gene_list:
		for i in runFirst:
			if i in gene:
				priority_gene_list.append(gene)
				remove.append(gene)
	for i in remove:
		gene_list.remove(gene)

#
#Running AnGST on unfinished gene trees
#
if priority_gene_list:
	p=Pool(cores)
	for gene in priority_gene_list:
		p.apply_async(runner, args=(gene,))
		pass
	p.close()
	p.join()

p=Pool(cores)
for gene in gene_list:
	p.apply_async(runner, args=(gene,))
	pass
p.close()
p.join()

#Compiling data from AnGST run
angst_data = os.listdir(angst_folder)
os.system("mkdir -p {0}/data_files".format(angst_folder))
out_sum = open("{0}/00_data_summary.out".format(angst_folder),'a')
summary = ["hgt inferred:", "dup inferred:", "los inferred:", "spc inferred:", "gene count:",
"run-time", "last reconciliation"]

for fi in angst_data:
	if "_angst" in fi:
		gene_data = os.listdir("{0}/{1}".format(angst_folder,fi))
		gene = fi
		gene = gene.replace("_angst", "")
		out = open("{0}/data_files/{1}_data".format(angst_folder,gene), 'w')
		out_sum.write("####    {0}    ####\n".format(gene))

		for f in gene_data:
			stats_flag = 0
			name = str(f)
			name = name.replace("AnGST.", "")
			if "stats" in name:
				stats_flag = 1

			out.write("\n#########\t{0}\t#########\n".format(name))
			data = open("{0}/{1}/{2}".format(angst_folder,fi,f),'r')

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

os.system("rm -f -r {0}/*_angst".format(angst_folder))
