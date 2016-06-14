from Bio import Phylo
from cStringIO import StringIO
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

#print alignment

fpin = open("scores.json",'r')
score_matrix = json.load(fpin)

count = 0
#clade_states = dict(zip(terms, [set([alignment[c.name]]) for c in terms])) # make dict of leaves and a set for TM
print clade_states 
for clade in tree.get_nonterminals(order="postorder"): #Gets internal nodes in a list going up the tree
    sub_branches = clade.clades #the children of each internal node
    #state1 = clade_states[sub_branches[0]] # get letter of each child
    #state2 = clade_states[sub_branches[1]]
    #state = state1 & state2 
    state1 = sub_branches[0].name
    state2 = sub_branches[1].name
    if score_matrix[state1+'\t'+state2] == 0 or score_matrix[state2+'\t'+state1] == 0: 
    	count += 1
    '''if not state: #if states not same, make set of two states and increment count
        state = state1 | state2
        count = count + 1
    clade_states[clade] = state
    print "internal node =", ' '.join(state)''' # print states and count so you can see how it traverses the tree
print "score =", count # count is equal to the number of mutations that have happened 



fpin.close()
