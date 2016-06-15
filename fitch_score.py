from Bio import Phylo
#from cStringIO import StringIO
import json

tree = Phylo.read('tree.txt', 'newick') #read from file
#tree = Phylo.read(StringIO("((A,G),T),(((G,G),A),((T,T),A));"), "newick") #Hardcoded as string

#if not tree.rooted:
	#tree.root_at_midpoint() #couldn't get this to work properly, supposed to root tree if unrooted

# sort tree terminals/leaves 
terms = tree.get_terminals()
terms.sort(key=lambda term: term.name)

# Dictionary of proteins and number of TMs
#alignment = {n.split(',')[0]: n.split(',')[1] for n in open("TMs.txt",'r').read().splitlines()}

fpin = open("scores.json",'r')
score_matrix = json.load(fpin) # load scores from json file 

count = 0
#clade_states = dict(zip(terms, [set([alignment[c.name]]) for c in terms])) # make dict of leaves and a set for TM

for clade in tree.get_nonterminals(order="postorder"): #Gets internal nodes in a list going up the tree
	sub_branches = clade.clades #the children of each internal node
	#state1 = clade_states[sub_branches[0]] # get letter of each child
	#state2 = clade_states[sub_branches[1]]
	#state = state1 & state2 
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
	"""iterate through each item the child nodes and check if combination is in dict and get value
	if there are no matching pairs, then iterate the count and label the next node with those names"""
	for n in state1: 
		for i in state2:
			if n+'\t'+i in score_matrix and score_matrix[n+'\t'+i] == 1:
				clade.name = [n]
				print 'match'
				match +=1
				break

			elif i+'\t'+n in score_matrix and score_matrix[i+'\t'+n] == 1:
				clade.name = [n]
				print 'match'
				match += 1
				break
					
	if match == 0:
		clade.name.extend(state1)
		clade.name.extend(state2)
		count+=1
		print "no match"

	print count 
	
	# from previous version that used number of TM segments 
	'''if not state: #if states not same, make set of two states and increment count
		state = state1 | state2
		count = count + 1
	clade_states[clade] = state
	print "internal node =", ' '.join(state)''' # print states and count so you can see how it traverses the tree

print "score =", count # count is equal to the number of mutations that have happened 



fpin.close()
