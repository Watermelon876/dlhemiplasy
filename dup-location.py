## get an output file with all the possible instances of hemiplasy in a data set

import os
import collections

from rasmus import util, treelib
from compbio import phylo

def hemiplasyConditions(numFamilies = 5351,
                        dataPath = '/home/muddcs15/research/work/hemiplasy/data/real-fungi/',
                        outputFile = '/home/muddcs15/research/work/hemiplasy/results/hemiplasy-loss.txt',
                        spectree = '/home/muddcs15/research/work/hemiplasy/data/config/fungi.stree'):
    # create variables and output file
    output = open(outputFile,'w')
    count = 0

    # define a list of all species and lists of each of the species pairs in separate lists
    stree = treelib.read_tree(spectree) # species tree
    species = stree.leaf_names()
    species1 = []
    species2 = []
    for node in stree:
        if len(node.leaves()) == 2:
            species1.append(node.children[0].name)
            species2.append(node.children[1].name)
    
    # loop over each fam id
    for famid in xrange(numFamilies):
        flag = False                                # this families met the criteria for possible hemiplasy
        locus_dict = collections.defaultdict(list)  # key = locus number, val = list of (gn, sp) in the locus
        famFilename = dataPath + '%d/%d-dup.dlcoal.dlcpar.recon' % (famid,famid)

        # if the file is not empty, process it
        if os.stat(famFilename).st_size != 0:
        
            # read the locus tree and the reconcilitation file
            tree_filename = dataPath + '%d/%d.dlcoal.locus.tree' % (famid,famid)
            recon_filename = dataPath + '%d/%d.dlcoal.locus.recon' % (famid,famid)
            tree = treelib.read_tree(tree_filename) # locus tree
            recon, events = phylo.read_recon_events(recon_filename, tree, stree) # reconciliation and events

            # create a dictionary for [locus] = species tree location
            locus_sname = {}
            
            # find location in species tree where each locus was created and then close dlcpar file
            dupFilename = dataPath + '%d/%d.dlcoal.dlcpar.dup.rel.txt' % (famid,famid)
            for line in util.open_stream(dupFilename):
                locus, gns1, gns2, sname = line.rstrip().split('\t')
                locus_sname[locus] = sname

            # track to genes and species in each locus
            for line in util.open_stream(famFilename):

                # assign names to the columns in the file
                gn, sp, locus = line.rstrip().split('\t')
                if locus == "1":
                    continue
                
                # store dict of key = locus, val = list of (gene, species) in locus
                locus_dict[locus].append((gn, sp))

            # for each locus, determine if genes in locus satisfy the properties for a possible hemiplasy
            for locus, lst in locus_dict.iteritems():
                sps = [sp for (gn,sp) in lst]
                gns = [gn for (gn,sp) in lst]
                
                # check if exists in only one species in a pair
                for sp1, sp2 in zip(species1, species2):
                    if (sp1 in sps and sp2 not in sps) or \
                        (sp2 in sps and sp1 not in sps):
                        
                        # check if exists elsewhere (outside pair)
                        for allsp in species:
                            if (allsp != sp1) and (allsp != sp2) and (allsp in sps):
                                flag = True

                # output this family id, the locus, the species with that locus, the genes on that locus,\
                # the species tree branch on which the duplication occurred, and the daughter of the duplication node in the locus tree
                if flag:
                    leaf_sps = []
                    leaf_gns = []
                    for gn, sp in lst:
                       if not gn.isdigit():
                           leaf_sps.append(sp)
                           leaf_gns.append(gn)

                    gnodes = [tree.nodes[name] for name in leaf_gns]
                    lca = treelib.lca(gnodes)
                    
                    output.write('\t'.join([str(famid), locus, ','.join(leaf_sps), ','.join(leaf_gns), locus_sname[locus], lca.name]))
                    output.write('\n')
                    break
            
            # if it is a true case, add to count
            if flag:
                count += 1
            
    # print total count and close output file        
    print "Total number of true cases =", count
    output.close()
    
if __name__ == "__main__":
    hemiplasyConditions()


    
