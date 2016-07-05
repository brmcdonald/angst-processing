import sys,os,re,math
from Bio import Phylo
import parser_functions as pf
from scipy import stats
import numpy as np

usage = """\n<Angst ouptut folder>\n<(D)omain or (P)rotein mode>\n<optional node naming file>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 2:
	angst_folder = argv.pop(0)
	mode_flag = argv.pop(0)
elif len(sys.argv) == 3:
	angst_folder = argv.pop(0)
	mode_flag = argv.pop(0)
	node_name_file = argv.pop(0)
else:
	sys.exit(usage)

#RE precompiles
get_tag_std = re.compile("([^\|]+)\|")
get_tag_locus = re.compile("([^_]+)_")
get_count = re.compile(":\s+(\d+)\n")

#Checks of NODE_N has an entry for a node yet, if not it makes one
def node_check(node):
	key = ''
	good_key = ''
	if not NODE_LOOKUP.has_key(node):
		term = re.split("-", node)
		if len(term) == 1:
			NODE_LOOKUP[term[0]] = term[0]
			good_key = term[0]
			NODE_N[good_key].angst_name = good_key
		else:
			found = 0
			for key in NODE_N:
				if len(NODE_N[key].term_list) == len(term):
					bad = 0
					for i in term:
						if not NODE_N[key].term_list.has_key(i):
							bad = 1
					if bad == 0:
						good_key = key
						NODE_N[good_key].angst_name = node
						NODE_LOOKUP[node] = good_key
						found = 1
			if found == 0:
				print "Cant find node for \n{0}".format(node)
				sys.exit()
	else:
		good_key = NODE_LOOKUP[node]
	if good_key:
		return good_key
	else:
		print 'Couldnt find good key for\n{0}'.format(node)
		sys.exit()

####################################################################################################
#######################################   Data Classes   ###########################################
####################################################################################################

#Container for data on a single gene family
class Fam():
	def __init__(self, name,lgt,dup,los,spec,total):
		self.name = name				#category name
		self.birth = ""					#name of birth node
		self.total = total				#total num genes
		self.lgt = lgt					#num of lgt events
		self.dup = dup					#num of dup events
		self.los = los					#num loss events
		self.spc = spc					#num speciation events
		self.event = lgt+dup+los+spc+1	#Total events
		self.g_contain = {}				#dict of genomes containing gene
		self.node_contain = {}			#dict of nodes containing gene
		self.genes = {}					#Genes in the family

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

#Contains data for each node of the species tree
class Node():
	def __init__(self, obj, flag):
		self.term_list = {}		#List of terminal nodes ie genomes
		self.dir_dec = {}		#List of direct decendent nodes
		self.dec = {}			#List of all decendent nodes
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

		#If its a genome, name it after the genome tag
		if self.genome == 1:
			name = self.node_obj.name
			m = get_tag_std.search(name)
			if m:
				self.name = m.group(1)
				self.term_list[m.group(1)] = ""
			else:
				self.name = name
				self.term_list[name] = ""

		#If not, give it a clade name and figure out terminal/decendent nodes
		if self.genome == 0:
			self.name = self.node_obj.name
			#Determines angst name using terminal nodes
			terminals = self.node_obj.get_terminals()
			for i in terminals:
				name = i.name
				m = get_tag_std.search(name)
				if m:
					self.term_list[m.group(1)] = ""
				else:
					self.term_list[name] = ""

			#Determines direct decendent nodes: get decendents and check them against each other
			subclades = list(self.node_obj.find_clades())
			nested = []
			for cladeA in subclades:
				if cladeA is self.node_obj:
					nested.append(cladeA)
				else:
					self.dec[cladeA.name] = ""		#Builds list of all decendent nodes
					for cladeB in subclades:
						if cladeA is not cladeB:
							if cladeA.is_parent_of(cladeB) and cladeA is not self.node_obj:
								nested.append(cladeB)
			for i in set(nested):
				subclades.remove(i)
			for i in subclades:
				self.dir_dec[i.name] = ""

class Gene():
	def __init__(self):
		self.name = ''
		self.fam = ''
		self.brn = ''
		self.genome = ''
		self.fam = ''
		self.dup = []
		self.lgt = []	#Each event is a list: [source,destination]
		self.spc = []
		self.event = []	#[source,type] of each event

####################################################################################################
###################################   Data Compilation   ###########################################
####################################################################################################


NODE_N = {}	#Node objects listed by node name
LEAF = {} 	#Keeps track of events that occur in the ancestors of each existing gene

#
#Data Summary parsing
#
data_sum = open("{0}/00_data_summary.out".format(angst_folder), 'r')
data_sum = data_sum.readlines()

#Compile basic statistics on each gene family
FAMS = {}	#Dict of category objects
re_count = re.compile("gene count:\s+(\d+)\n")
re_lgt = re.compile("hgt inferred:\s+(\d+)\n")
re_dup = re.compile("dup inferred:\s+(\d+)\n")
re_los = re.compile("los inferred:\s+(\d+)\n")
re_spc = re.compile("spc inferred:\s+(\d+)\n")
data_res = [re_count,re_lgt,re_dup,re_los,re_spc]

total_lgt = 0
total_dup = 0
total_los = 0
total_spc = 0
total_gene_current = 0
total_gene_treewide = 0
total_event = 0
for ln in data_sum:
	if "####    " in ln:
		ln = ln.replace("####    ", "")
		name = ln.replace("    ####\n", "")
		name = re.sub("_boot\S+\d+$","",name)
		lgt = 0
		count = 0
		dup = 0
		los = 0
		spc = 0

	elif "run-time" in ln:
		FAMS[name] = Fam(name,lgt,dup,los,spc,count)
		name = ""
	else:
		for idx, data_re in enumerate(data_res):
			m = data_re.search(ln)
			if m and idx == 0:
				count = float(m.group(1))
				total_gene_current += float(m.group(1))
				total_event += float(m.group(1))
			elif m and idx == 1:
				lgt = float(m.group(1))
				total_lgt += float(m.group(1))
				total_event += float(m.group(1))
			elif m and idx == 2:
				dup = float(m.group(1))
				total_dup += float(m.group(1))
				total_event += float(m.group(1))
			elif m and idx == 3:
				los = float(m.group(1))
				total_los += float(m.group(1))
				total_event += float(m.group(1))
			elif m and idx == 4:
				spc = float(m.group(1))
				total_spc += float(m.group(1))
				total_event += float(m.group(1))

#
#Species tree and node parsing
#
species_tree = open("{0}/species_tree_clean".format(angst_folder), 'r')
species_tree = species_tree.readlines()
out_temp = open("{0}/species_tree.newick".format(angst_folder), 'w')
for ln in species_tree:
	ln = ln.replace("/", "_")
	ln = list(ln)
	ln.pop()
	ln.pop()
	ln.pop(0)
	ln = ''.join(ln)
	out_temp.write(ln)
	out_temp.write(";")
out_temp.close()

#Loading tree into BioPython and naming any unnamed nodes
label = 0
inTree = open("{0}/species_tree.newick".format(angst_folder), 'r')
inTree = inTree.readlines()
output = open("{0}/species_tree_labels.newick".format(angst_folder), "w")
for ln in inTree:
	ln = ln.translate(None, "\n\r")
	ln = ln.split("):")
	tail = ln.pop()
	for i in ln:
		output.write("{0})NODE{1}:".format(i,label))
		label += 1
	output.write(tail)
output.close()

tree = Phylo.read("{0}/species_tree_labels.newick".format(angst_folder),'newick')
genome_lst = tree.get_terminals()
node_lst = tree.get_nonterminals()

total_nodes = len(genome_lst) + len(node_lst)

#
#Build tree structure: each node has a list of direct decendents and a parent node
#

#Builds a dictionary of node objects based on node name.
for idx, i in enumerate(node_lst):
	node = Node(i,0)
	NODE_N[node.name] = node
for i in genome_lst:
	node = Node(i,1)
	NODE_N[node.name] = node

#Identifies parent node for each node
for i in NODE_N:
	dec = NODE_N[i].dir_dec
	for j in dec:
		m = get_tag_std.search(j)
		if m:
			NODE_N[m.group(1)].parent = i
		else:
			NODE_N[j].parent = i

#
#Node data parsing
#
reNODE = re.compile("\:\s+(.+)\n")
gene_files = os.listdir("{0}/data_files".format(angst_folder))

NODE_LOOKUP = {} #lookup table: AnGST node name -> Node name

#Initializing dicts for each node
for f in gene_files:
	g_name = str(f)
	g_name = re.sub("_boot\S+", "", g_name)

	for node in NODE_N:
		NODE_N[node].dup[g_name] = 0
		NODE_N[node].los[g_name] = 0
		NODE_N[node].spc[g_name] = 0
		NODE_N[node].cats[g_name] = 0
		NODE_N[node].lgt_target_gene[g_name] = 0
		NODE_N[node].lgt_source_gene[g_name] = 0
		NODE_N[node].event[g_name] = 0

#Collects data on events at each node
if mode_flag == 'P':
	for f in gene_files:
		g_name = str(f)
		g_name = re.sub("_boot\S+", "", g_name)

		with open("{0}/data_files/{1}".format(angst_folder,f),'r') as f_in:
			data_type = None
			for ln in f_in:
				if "#########	counts	#########" in ln:
					data_type = "counts"
				elif "#########	events	#########" in ln:
					data_type = "events"
				elif '#########	leaf	#########' in ln:
					data_type = 'leaf'
					x = ''
				elif "############################" in ln:
					if data_type == 'leaf':
						LEAF[x.genome][x.name] = x
					data_type = None

				elif ":" in ln and data_type:
					if data_type == "counts":
						m = re.search("\:\s+(\d+)\n$", ln)
						if m:
							ln_clean = ln.replace(m.group(0), "")
							count = m.group(1)

							#If the angst node name has already been assigned, sort data
							if NODE_N.has_key(ln_clean):
								NODE_N[ln_clean].total += int(count)
								total_gene_treewide += int(count)
								NODE_N[ln_clean].cats[g_name] += int(count)

							#If not, identify the correct angst node name by finding terminal nodes
							else:
								node_check(ln_clean)
								ln_clean = NODE_LOOKUP[ln_clean]
								NODE_N[ln_clean].total += int(count)
								total_gene_treewide += int(count)
								NODE_N[ln_clean].cats[g_name] += int(count)

					#Data tracing the history of each extant gene
					elif data_type == 'leaf':
						if 'leaf ' in ln:
							if x:
								LEAF[x.genome][x.name] = x
							ln = ln.replace(':\n','')
							ln = ln.split(' ')
							x = Gene()
							x.name = ln[1]
							x_name = ln[1].split('_')
							x.genome = x_name[0]
							x.fam = g_name
							FAMS[g_name].genes[ln[1]] = ''
							if not LEAF.has_key(x.genome):
								LEAF[x.genome] = {}
						elif '[dup]' in ln:
							ln = ln.replace('\n','')
							ln = ln.split(': ')
							if '-' in ln[1]:
								node = node_check(ln[1])
							else:
								node = ln[1]
							x.dup.append(node)
							x.event.append(['dup',node])
						elif '[hgt]' in ln:
							ln = ln.replace('\n','')
							ln = ln.split(': ')
							data = ln[1].split(' --> ')
							if '-' in data[0]:
								node_a = node_check(data[0])
							else:
								node_a = data[0]
							if '-' in data[1]:
								node_b = node_check(data[1])
							else:
								node_b = data[1]
							x.lgt.append([node_a,node_b])
							x.event.append(['lgt',node_a])
						elif '[spc]' in ln:
							ln = ln.replace('\n','')
							ln = ln.split(': ')
							if '-' in ln[1]:
								node = node_check(ln[1])
							x.spc.append(node)
							x.event.append(['spc',node])

					#Events data:births,deaths,lgt,duplications
					elif data_type == "events":
						m = reNODE.search(ln)

						if '[hgt]' in ln:
							node_1 = ''
							node_2 = ''
							if m:
								node = re.split("\s+\-\-\>\s+", m.group(1))
								node_1 = node_check(node[0])
								node_2 = node_check(node[1])

								#Record source node for lgt
								NODE_N[node_1].lgt_total_source += 1
								NODE_N[node_1].lgt_source_gene[g_name] += 1
								if NODE_N[node_1].lgt_source.has_key(node_2):
									NODE_N[node_1].lgt_source[node_2].append(g_name)
								else:
									NODE_N[node_1].lgt_source[node_2] = []
									NODE_N[node_1].lgt_source[node_2].append(g_name)

								#All other data relates to the recieving node
								NODE_N[node_2].lgt_total_target += 1
								NODE_N[node_2].lgt_target_gene[g_name] += 1
								NODE_N[node_2].event_total += 1
								NODE_N[node_2].event[g_name] += 1

								if NODE_N[node_2].lgt_target.has_key(node_1):
									NODE_N[node_2].lgt_target[node_1].append(g_name)
								else:
									NODE_N[node_2].lgt_target[node_1] = []
									NODE_N[node_2].lgt_target[node_1].append(g_name)
							else:
								print 'Cant read hgt data for {0}'.format(ln)

						elif "[dup]" in ln:
							node = node_check(m.group(1))
							NODE_N[node].dup[g_name] += 1
							NODE_N[node].dup_total += 1
							NODE_N[node].event_total += 1
							NODE_N[node].event[g_name] += 1
						elif "[los]" in ln:
							node = node_check(m.group(1))
							NODE_N[node].los[g_name] += 1
							NODE_N[node].los_total += 1
							NODE_N[node].event_total += 1
							NODE_N[node].event[g_name] += 1
						elif "[brn]" in ln:
							node = node_check(m.group(1))
							NODE_N[node].brn_total += 1
							NODE_N[node].brn[g_name] = 1
							NODE_N[node].event_total += 1
							NODE_N[node].event[g_name] += 1

#For parsing domain gene names
elif mode_flag == 'D':
	for f in gene_files:
		g_name = str(f)
		g_name = re.sub("_boot\S+", "", g_name)

		with open("{0}/data_files/{1}".format(angst_folder,f),'r') as f_in:
			data_type = None
			for ln in f_in:
				if "#########	counts	#########" in ln:
					data_type = "counts"
				elif "#########	events	#########" in ln:
					data_type = "events"
				elif '#########	leaf	#########' in ln:
					data_type = 'leaf'
					x = ''
				elif "############################" in ln:
					if data_type == 'leaf':
						LEAF[x.genome][x.name] = x
					data_type = None

				elif ":" in ln and data_type:
					if data_type == "counts":
						m = re.search("\:\s+(\d+)\n$", ln)
						if m:
							ln_clean = ln.replace(m.group(0), "")
							count = m.group(1)

							#If the angst node name has already been assigned, sort data
							if NODE_N.has_key(ln_clean):
								NODE_N[ln_clean].total += int(count)
								total_gene_treewide += int(count)
								NODE_N[ln_clean].cats[g_name] += int(count)

							#If not, identify the correct angst node name by finding terminal nodes
							else:
								node_check(ln_clean)
								ln_clean = NODE_LOOKUP[ln_clean]
								NODE_N[ln_clean].total += int(count)
								total_gene_treewide += int(count)
								NODE_N[ln_clean].cats[g_name] += int(count)

					elif data_type == 'leaf':
						if 'leaf ' in ln:
							if x:
								LEAF[x.genome][x.name] = x
							ln = ln.replace(':\n','')
							ln = ln.split(' ')
							x = Gene()
							x.name = ln[1]
							x_name = ln[1].split('_')
							x.genome = x_name[0]
							x.fam = g_name
							if not LEAF.has_key(x.genome):
								LEAF[x.genome] = {}
						elif '[dup]' in ln:
							ln = ln.replace('\n','')
							ln = ln.split(': ')
							if '-' in ln[1]:
								node = node_check(ln[1])
							else:
								node = ln[1]
							x.dup.append(node)
						elif '[hgt]' in ln:
							ln = ln.replace('\n','')
							ln = ln.split(': ')
							data = ln[1].split(' --> ')
							if '-' in data[0]:
								node_a = node_check(data[0])
							else:
								node_a = data[0]
							if '-' in data[1]:
								node_b = node_check(data[1])
							else:
								node_b = data[1]
							x.lgt.append([node_a,node_b])
						elif '[spc]' in ln:
							ln = ln.replace('\n','')
							ln = ln.split(': ')
							if '-' in ln[1]:
								node = node_check(ln[1])
							x.spc.append(node)

					#Events data:births,deaths,lgt,duplications
					elif data_type == "events":
						m = reNODE.search(ln)

						if '[hgt]' in ln:
							node_1 = ''
							node_2 = ''
							if m:
								node = re.split("\s+\-\-\>\s+", m.group(1))
								node_1 = node_check(node[0])
								node_2 = node_check(node[1])

								#Record source node for lgt
								NODE_N[node_1].lgt_total_source += 1
								NODE_N[node_1].lgt_source_gene[g_name] += 1
								if NODE_N[node_1].lgt_source.has_key(node_2):
									NODE_N[node_1].lgt_source[node_2].append(g_name)
								else:
									NODE_N[node_1].lgt_source[node_2] = []
									NODE_N[node_1].lgt_source[node_2].append(g_name)

								#All other data relates to the recieving node
								NODE_N[node_2].lgt_total_target += 1
								NODE_N[node_2].lgt_target_gene[g_name] += 1
								NODE_N[node_2].event_total += 1
								NODE_N[node_2].event[g_name] += 1

								if NODE_N[node_2].lgt_target.has_key(node_1):
									NODE_N[node_2].lgt_target[node_1].append(g_name)
								else:
									NODE_N[node_2].lgt_target[node_1] = []
									NODE_N[node_2].lgt_target[node_1].append(g_name)
							else:
								print 'Cant read hgt data for {0}'.format(ln)

						elif "[dup]" in ln:
							node = node_check(m.group(1))
							NODE_N[node].dup[g_name] += 1
							NODE_N[node].dup_total += 1
							NODE_N[node].event_total += 1
							NODE_N[node].event[g_name] += 1
						elif "[los]" in ln:
							node = node_check(m.group(1))
							NODE_N[node].los[g_name] += 1
							NODE_N[node].los_total += 1
							NODE_N[node].event_total += 1
							NODE_N[node].event[g_name] += 1
						elif "[brn]" in ln:
							node = node_check(m.group(1))
							NODE_N[node].brn_total += 1
							NODE_N[node].brn[g_name] = 1
							NODE_N[node].event_total += 1
							NODE_N[node].event[g_name] += 1

#Add 'spc' events for existing genes to genome nodes, which dont get any spc events
for i in NODE_N:
	NODE_N[i].spc_total = NODE_N[i].total - (NODE_N[i].lgt_total_target + (NODE_N[i].dup_total)) - NODE_N[i].brn_total
	NODE_N[i].event_total = NODE_N[i].event_total + NODE_N[i].spc_total
	total_spc += NODE_N[i].spc_total
	total_event += NODE_N[i].spc_total
	for fam in NODE_N[i].cats:
		adjust = abs(NODE_N[i].lgt_target_gene[fam] + (NODE_N[i].dup[fam]))

		NODE_N[i].spc[fam] = NODE_N[i].cats[fam] - adjust
		NODE_N[i].event[fam] = NODE_N[i].cats[fam] - NODE_N[i].event[fam]

#Importing global data to parsing script
pf.import_data(LEAF,NODE_N,FAMS,tree)

####################################################################################################
#######################################   Analysis   ###############################################
####################################################################################################



#lgt_tree("./event_lgt.itol")
#pf.piechart_tree_proportion(NODE_N,"./event_piechart.itol")
#pf.cat_data_table(FAMS,"./family_event_summary.out", "N")

#DATA = pf.build_kegg("k2c")
#CATS = DATA[1]
#KO = DATA[0]

#pf.fet_cat(CATS, NODE_N,total_lgt,total_gene_treewide,"lgt")
#pf.fet_fam("lgt")

#ISO = pf.node_isolation(NODE_N,"strept_isolation","itol",'NODE171',tree)
#event_rate('brn','itol')

