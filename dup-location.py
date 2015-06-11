## get an output file with all the possible instances of hemiplasy in a data set

import os
import collections

from rasmus import util, treelib
from compbio import phylo

stree = treelib.read_tree('/home/muddcs15/research/work/hemiplasy/data/config/fungi.stree') # species tree

def hemiplasy():
    # create variables and output file
    output = open('/home/muddcs15/research/work/hemiplasy/data/real-fungi-rel/hemiplasy-loss.txt','w')
    count = 0

    # define a list of all species and lists of each of the species pairs in separate lists
    species = ['scer','agos','calb','cpar','cgui','spar','klac''ctro','lelo','dlhan','smik','sbay','cgla','scas','kwal','clus']
    species1 = ['scer','agos','calb','cpar','cgui']
    species2 = ['spar','klac','ctro','lelo','dlhan']
    
    # loop over each fam id
    for famid in xrange(5351):
        # make standard check = False, create locus dictionary, convert famid to a string, and open dup file of fam id
        check = False
        locus_dict = collections.defaultdict(list)
        strfam = str(famid)
        famFile = open('/home/muddcs15/research/work/hemiplasy/data/real-fungi/%d/%d-dup.dlcoal.dlcpar.recon' % (famid,famid))
        statinfo = os.stat('/home/muddcs15/research/work/hemiplasy/data/real-fungi/%d/%d-dup.dlcoal.dlcpar.recon' % (famid,famid))
        
        # read the locus tree and the reconcilitation file
        tree_filename = '/home/muddcs15/research/work/hemiplasy/data/real-fungi/%d/%d.dlcoal.locus.tree' % (famid,famid)
        recon_filename = '/home/muddcs15/research/work/hemiplasy/data/real-fungi/%d/%d.dlcoal.locus.recon' % (famid,famid)
        tree = treelib.read_tree(tree_filename) # locus tree
        recon, events = phylo.read_recon_events(recon_filename, tree, stree) # reconciliation and events

        # if the file is not empty, run the loop
        if statinfo.st_size != 0:

            # create a dictionary for [locus] = species tree location
            locus_sname = {}
            f = open('/home/muddcs15/research/work/hemiplasy/data/real-fungi/%d/%d.dlcoal.dlcpar.dup.rel.txt' % (famid,famid))
            # find location in species tree where each locus was created
            for line in f:
                locus, gns1, gns2, sname = line.rstrip().split('\t')
                locus_sname[locus] = sname
            f.close()

            # loop through each line in the dlcpar result of each family id
            for line in famFile:

                # assign names to the columns in the file
                gn, sp, locus = line.rstrip().split('\t')
                
                # store dict of key = locus, val = list of (gene, species) in locus
                locus_dict[locus].append((gn, sp))
                
                # for each locus, determine if genes in locus satisfy the properties for a possible hemiplasy
                sps = [sp for (gn,sp) in lst]
                gns = [gn for (gn,sp) in lst]
                
                # check if exists in only one species in a pair
                for sp1, sp2 in zip(species1, species2):
                    if (sp1 in sps and sp2 not in sps) or \
                        (sp2 in sps and sp1 not in sps):
                        
                        # check if exists elsewhere (outside pair)
                        for allsp in species:
                            if (allsp != sp1) and (allsp != sp2) and (allsp in sps):
                                check = True

                # output this family id, the locus, the species with that locus, the genes on that locus,\
                # the species tree branch on which the duplication occurred, and the daughter of the duplication node in the locus tree
                if check:
                    leaf_sps = []
                    leaf_gns = []
                    for gn, sp in lst:
                       if not gn.isdigit():
                           leaf_sps.append(sp)
                           leaf_gns.append(gn)

                    gnodes = [tree.nodes[name] for name in leaf_gns]
                    lca = treelib.lca(gnodes)
                    
                    output.write('\t'.join([strfam, locus, ','.join(leaf_sps), ','.join(leaf_gns), locus_sname[locus], lca.name]))
                    output.write('\n')
                    break

        # close file
        famFile.close()
        
        # if it is a true case, add to count
        if check:
            count += 1
            
    # print total count and close output file        
    print "Total number of true cases =", count
    output.close()
    
if __name__ == "__main__":
    hemiplasy()

    


    
