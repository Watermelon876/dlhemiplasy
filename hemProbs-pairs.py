import os
import matplotlib.pyplot as plt
from rasmus import stats, treelib

def boxPlot(dataPath = '/home/muddcs15/research/work/hemiplasy/results/',
            prob1 = '0.001',
            prob2 = '0.05',
            prob3 = '0.1',
            prob4 = '0.5',
            spectree = '/home/muddcs15/research/work/hemiplasy/data/config/fungi.stree'):
    """
    A function that will output boxplots of probability of hemiplasy and probability of hemiplasy over duploss vs. initial allele frequency
    """
    

    stree = treelib.read_tree(spectree) # species tree
    species = stree.leaf_names()
    species1 = []
    species2 = []
    for node in stree:
        if len(node.leaves()) == 2:
            species1.append(node.children[0].name)
            species2.append(node.children[1].name)

    # identify the files for each of the different initial frequencies
    probs1 = os.path.join(dataPath, 'probabilities-' + prob1 + '.txt')
    probs2 = os.path.join(dataPath, 'probabilities-' + prob2 + '.txt')
    probs3 = os.path.join(dataPath, 'probabilities-' + prob3 + '.txt')
    probs4 = os.path.join(dataPath, 'probabilities-' + prob4 + '.txt')
    probsList = [probs1, probs2, probs3, probs4]

    totalPerList = []   # probability of hemiplasy compared to duploss
    totalAList = []     # probability of hemiplasy ocurring
    totalPairs = []

    h = 0 # probability that ocurred by hemiplasy
    d = 0 # probability that ocurred by duploss

    pair1 = []
    pair2 = []
    pair3 = []
    pair4 = []
    pair5 = []

    # open each probability file
    for probFilename in probsList:
        events = open('/home/muddcs15/research/work/hemiplasy/results/hemiplasy-loss.txt', 'r')
        
        hList = []      # list of probability of hemiplasy
        perList = []    # list of percentage with prob hemiplasy > prob duploss
        aveList = []    # list of average probability of hemiplasy per fam id

        countTrue = 0

        # look at each famid for that initial frequency
        probFile = open(probFilename, "r")
        for line in probFile:
            sepProbs = line.split()
            fam = sepProbs.pop(0)
            famid = fam[6:]

            # get the probability of duploss and hemiplasy for each trial in each famid
            for pair in sepProbs:
                duploss, hemiplasy = map(float, pair.split(','))
                hList.append(hemiplasy)
                # check whether hemiplasy is more likely or duploss
                if hemiplasy > duploss:
                    h += 1
                else:
                    d += 1

            # calculate the percent that likely ocurred by hemiplasy            
            percent = float(h)/float(h+d)
            # get the average probability of hemiplasy for each famid
            ave = stats.mean(hList)

            # append percent by hemiplasy to perList and average for the famid to aveList
            perList.append(percent)
            aveList.append(ave)
            
            for line in events:
                ev_famid, locus, spcs, gns, snode, lca = line.rstrip().split('\t')
                if famid == ev_famid:
                    countTrue += 1
                    for sp1, sp2 in zip(species1, species2):
                        if (sp1 in spcs and sp2 not in spcs):
                            spec_check = sp1
                            specPos = species1.index(sp1)
                            
                        elif (sp2 in spcs and sp1 not in spcs):
                            spec_check = sp2
                            specPos = species2.index(sp2)
                    break
                    
            if specPos == 0:
                pair1.append(ave)
            if specPos == 1:
                pair2.append(ave)
            if specPos == 2:
                pair3.append(ave)
            if specPos == 3:
                pair4.append(ave)
            if specPos == 4:
                pair5.append(ave)
        events.close()
                
        # append the lists through each famid to the large lists for each list of values
        totalPerList.append(perList)
        
        totalAList.append(aveList)
        

        # close file
        probFile.close()
        
    totalPairs.append(pair1)
    totalPairs.append(pair2)
    totalPairs.append(pair3)
    totalPairs.append(pair4)
    totalPairs.append(pair5)
    
    plt.boxplot(totalPairs)
    plt.title('Hemiplasy by Pairs')
    plt.xlabel('Pair')
    plt.ylabel('Probability')
    
    # print the plots
    plt.show()

if __name__ == "__main__":
    boxPlot()
