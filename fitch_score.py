#!/usr/bin/env python

from Bio import Phylo
from Bio import SeqIO
import re

tree = Phylo.read('tree.txt', 'newick') #read from file

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

def GetTMPosition_gapless(topo):#{{{
    """
    Get the position of TM helices given the topology (without gaps)
    The return value is a list of 2-tuples: [ (beg, end), (beg, end)...]
    """
    posTM=[]
    m=re.finditer("(M+)",topo)
    for i in m:
        posTM.append((i.start(0), i.end(0)))
    return posTM

def CountTM(topo):#{{{
    """Count the number of TM regions in a topology with or without gaps"""
    return len(GetTMPosition_gapless(topo))

#make dictionary of protein IDs and their topologies from topology file using SeqIO
topos = {}
for seq_record in SeqIO.parse("TMs.txt", "fasta"):
	if GetNtermState(seq_record.seq) == 'i':
		nterm = '-'
	else:
		nterm = '+'
	topos[seq_record.id[3:9]] = nterm + str(CountTM(str(seq_record.seq)))
	
#print topos

count = 0

for clade in tree.get_nonterminals(order="postorder"): #Gets internal nodes in a list going up the tree
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
	print state1, state2

	clade.name = []
	match = 0
	"""iterate through each item the child nodes and check if their topologies are the same
	if there are no matching pairs, then iterate the count and label the next node with those names"""
	for n in state1: 
		for i in state2:
			print topos[n], topos[i]
			if topos[n] == topos[i]:
				clade.name = [n]
				print 'match'
				match +=1
				break
					
	if match == 0:
		clade.name.extend(state1)
		clade.name.extend(state2)
		count+=1
		print "no match"
	
	print "score =", count # count is equal to the number of mutations that have happened 


