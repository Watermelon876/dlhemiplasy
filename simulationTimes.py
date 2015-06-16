
from rasmus import util, treelib
from compbio import phylo
from dlcoal import duploss

import StringIO
import re

def times(birthRate = 0.000732,
          deathRate = 0.000859,
          premn = 90.0,
          trials = 100,
          dataPath = '/home/muddcs15/research/work/hemiplasy/data/'):
    """
    A funtion that outputs the time from a duplication to a node
    INPUT:
    birthrate
    deathrate  
    premn = float value of premean
    trials = number of times to run sample_dup_times to get multiple differences between speciation and duplication events
    """
    # open the file with the list of duplications and read the species tree
    duplications = open('/home/muddcs15/research/work/hemiplasy/results/hemiplasy-loss.txt','r')
    stree = treelib.read_tree(dataPath + 'config/fungi.stree')

    # get the ages of the nodes and the timestamps of the species tree
    times = treelib.get_tree_timestamps(stree)
    ages = treelib.get_tree_ages(stree)

    # open an output file
    output = open('/home/muddcs15/research/work/hemiplasy/input/totaltimes1.txt','w')

    # go through each family id with a duplication that met the prerequisites
    for line in duplications:
        
        # assign names to each of the columns in the duplications file
        famid, locus, leaf_sps, leaf_gns, dup_sp, daughter_snode_sp = line.rstrip().split('\t')
                
        # read the locus tree and the reconcilitation file
        locus_tree_filename = dataPath + 'real-fungi/%s/%s.dlcoal.locus.tree' % (famid, famid)
        locus_recon_filename = dataPath + 'real-fungi/%s/%s.dlcoal.locus.recon' % (famid, famid)
        ltree = treelib.read_tree(locus_tree_filename)
        recon, events = phylo.read_recon_events(locus_recon_filename, ltree, stree)

        # read daughter node from file
        daughter_node = ltree.nodes[daughter_snode_sp]
        node = daughter_node.parent

        # get age of the speciation node
        snode = stree.nodes[int(dup_sp)]
        nodeAge = ages[snode]
           
        # create an empty list to input times
        time_list=[]
        
        # run the number of inputted trials to get a list of time difference between speciation and duplication events
        for i in range(trials):
            nodeTime = duploss.sample_dup_times(ltree, stree, recon, birthRate, deathRate, pretime=None, premean=premn)
            t_x=nodeTime[node]
            dupAge = t_x - float(nodeAge)
            # make sure there is no computational error when calculating dupAge
            if (0 <= dupAge < 1e-13) or (0 >= dupAge > -1e-13) :
                dupAge = 0.0
            # add the dupAge to the time list     
            time_list+=[str(dupAge)]
            
        # convert the time list to comma separated string
        time_string = ",".join(time_list)

        # output the family id, the node age, and the difference between duplication time and node to output file
        output.write(famid + '\t' + str(nodeAge) + '\t' + time_string + '\n')

    # close the opened files
    duplications.close()
    output.close()
    
    
if __name__ == "__main__":
    times()
