#!/usr/bin/env python

from Bio import Phylo
from Bio import SeqIO
import re
import sys


usage = """
usage: %s treefile topofile
"""%(sys.argv[0])

if len(sys.argv) != 3:
    print usage
    sys.exit(1)

try:
    treefile = sys.argv[1]
    topofile = sys.argv[2]
except IndexError:
    print usage
    sys.exit(1)

    

# topology related functions from Nanjiang
GAP = '-'
def GetNtermState(topo):#{{{
    if topo[0] != GAP:
        return topo[0]
    else:
        topo=topo.lstrip(GAP)
        if topo != "":
            return topo[0]
        else:
            return None
 #}}} 

def GetTMPosition(topo):#{{{
    """
    Get position of TM helices given a topology
    this version is much faster (~25 times) than using than finditer
    updated 2011-10-24
    """
    posTM=[]
    lengthTopo=len(topo)
    b=0
    e=0
    while 1:
        b=topo.find('M',e)
        if b != -1:
            m = re.search('[io]', topo[b+1:])
            if m != None:
                e = m.start(0)+b+1
            else:
                e=lengthTopo
            if topo[e-1] == GAP:
                e=topo[:e-1].rfind('M')+1
            if b == e:
                print "Error topo[b-10:e+10]=", topo[b-30:e+30]
                #sys.exit(1)
                return []
            posTM.append((b,e))
        else:
            break
    return posTM
 #}}} 

def CountTM(topo):#{{{
    """Count the number of TM regions in a topology with or without gaps"""
    return len(GetTMPosition(topo))
 #}}} 


def FitchScore_snow1(treefile, topofile):#{{{
	"""
	Snow's old implementation
	"""
	tree = Phylo.read(treefile, 'newick') #read from file into tree object
#make dictionary of protein IDs and their topologies from topology file using SeqIO
	topos = {}
	for seq_record in SeqIO.parse(topofile, "fasta"):
		if GetNtermState(seq_record.seq) == 'i':
			nterm = '-'
		else:
			nterm = '+'
		topos[seq_record.id[3:9]] = nterm + str(CountTM(str(seq_record.seq)))

	count = 0
	nodes = tree.get_nonterminals(order="postorder")#Gets internal nodes in a list going up the tree
	for clade in nodes: 
		sub_branches = clade.clades #the children of each internal node
		print "sub_branches=",sub_branches
		state1 = [] # I want these as lists so I can iterate through them 
		state2 = []
		print "state1=", state1
		print "state2=", state2
		if type(sub_branches[0].name) == str: # This block is pretty ugly but I couldn't find a 
			state1.append(sub_branches[0].name) # better way of dealing with the single strings at the beginning
		if type(sub_branches[1].name) == str:
			state2.append(sub_branches[1].name)
		if type(sub_branches[0].name) == list:
			state1 = sub_branches[0].name
		if type(sub_branches[1].name) == list:
			state2 = sub_branches[1].name

		clade.name = []
		match = 0
		"""iterate through each item the child nodes and check if their topologies are the same;
		if there are no matching pairs, then iterate the count and label the next node with those names"""
		for n in state1: 
			for i in state2:
				if topos[n] == topos[i]:
					clade.name = [n]
					match +=1

		if match == 0:
# 		print "clade.name before extend=", clade.name
			clade.name.extend(state1)
			clade.name.extend(state2)
# 		print "clade.name after extend=", clade.name
			print
			count+=1
# 	else:
# 	    print "clade.name with match=", clade.name
		
	print "The total number of mutation events is:", count # count is equal to the number of mutations that have happened 

	'''Iterate throught the nodes in the tree and remove child nodes if they are same topology as parent so that 
	it leaves only the changes'''
	for clade in nodes:
		if type(clade.name) == list:
			clade.name = ','.join(clade.name) 
		for child in clade.clades:
			if child.name in topos and clade.name in topos:
				if topos[child.name] == topos[clade.name]:
					try:
						tree.prune(target=child)
					except ValueError:
						pass

		# This renames the internal nodes that are left with their topologies
		clade.name = ', '.join([topos[n] for n in clade.name.split(',')]) 

# This renames the terminals with their topologies 
	for term in tree.get_terminals():
		term.name = topos[term.name]
					
	tree.ladderize()
	Phylo.draw(tree, show_confidence=False) # can be changed to draw_ascii and also to show_confidence=False

#}}}

def FitchScore_nj1(treefile, topofile):#{{{
	"""
	Method by Nanjiang
	"""
	tree = Phylo.read(treefile, 'newick') #read from file into tree object
#make dictionary of protein IDs and their topologies from topology file using SeqIO
	topos = {}
	for seq_record in SeqIO.parse(topofile, "fasta"):
		if GetNtermState(seq_record.seq) == 'i':
			nterm = '-'
		else:
			nterm = '+'
		topos[seq_record.id[3:9]] = nterm + str(CountTM(str(seq_record.seq)))

	count = 0
	nodes = tree.get_nonterminals(order="postorder")#Gets internal nodes in a list going up the tree
	# first replace the name of all leafs by its state (nterm+numTM)

	# first, convert the label of each leaf to its state to be used for comparison
	for clade in nodes:
		sub_branches = clade.clades
		if type(sub_branches[0].name) == str: 
			sub_branches[0].name = topos[sub_branches[0].name]
		if type(sub_branches[1].name) == str:
			sub_branches[1].name = topos[sub_branches[1].name]
		if type(sub_branches[0].name) == list:
			sub_branches[0].name = [topos[x] for x in sub_branches[0].name]
		if type(sub_branches[1].name) == list:
			sub_branches[1].name = [topos[x] for x in sub_branches[1].name]

	for clade in nodes: 
		sub_branches = clade.clades #the children of each internal node
		state1 = [] # I want these as lists so I can iterate through them 
		state2 = []
		if type(sub_branches[0].name) == str: # This block is pretty ugly but I couldn't find a 
			state1.append(sub_branches[0].name) # better way of dealing with the single strings at the beginning
		if type(sub_branches[1].name) == str:
			state2.append(sub_branches[1].name)
		if type(sub_branches[0].name) == list:
			state1 = sub_branches[0].name
		if type(sub_branches[1].name) == list:
			state2 = sub_branches[1].name

		state_inter = set(state1) & set(state2)
		state_union = set(state1) | set(state2)

		match = len(state_inter)

		clade.name = []
		if match == 0:
			clade.name = list(state_union)
			count += 1
		else:
			clade.name = list(state_inter)

	print "The total number of mutation events is:", count # count is equal to the number of mutations that have happened 

	'''Iterate throught the nodes in the tree and remove child nodes if they are same topology as parent so that 
	it leaves only the changes'''
	for clade in nodes:
		if type(clade.name) == list:
			clade.name = ','.join(clade.name) 
		for child in clade.clades:
			if child.name == clade.name:
				try:
					tree.prune(target=child)
				except ValueError:
					pass

		# This renames the internal nodes that are left with their topologies
		clade.name = ', '.join(clade.name.split(','))

# This renames the terminals with their topologies 
					
	tree.ladderize()
	Phylo.draw(tree, show_confidence=False) # can be changed to draw_ascii and also to show_confidence=False#}}}



FitchScore_snow1(treefile, topofile)
FitchScore_nj1(treefile, topofile)

