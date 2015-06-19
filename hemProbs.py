import os
import collections
import matplotlib.pyplot as plt
from rasmus import stats, treelib

def boxPlot(x_axis = 'initialFreq',
            dataPath = '/home/muddcs15/research/work/hemiplasy/results/',
            prob1 = '0.001',
            prob2 = '0.01',
            prob3 = '0.05',
            prob4 = '0.1',
            prob5 = '0.5',
            spectree = '/home/muddcs15/research/work/hemiplasy/data/config/fungi.stree'):
    """
    A function that will output boxplots of probability of hemiplasy and probability of hemiplasy over duploss vs. initial allele frequency
    """

    if x_axis == 'initialFreq':
        events = open('/home/muddcs15/research/work/hemiplasy/results/hemiplasy-loss.error.txt', 'r')

        for line in events:
            famid, locus, spcs, gns, snode, lca = line.rstrip().split('\t')

        # define number of plots to be outputed
        fig, axes = plt.subplots(nrows=1, ncols=2)

        # identify the files for each of the different initial frequencies
        probs1 = os.path.join(dataPath, 'probabilities-' + prob1 + '.txt')
        probs2 = os.path.join(dataPath, 'probabilities-' + prob2 + '.txt')
        probs3 = os.path.join(dataPath, 'probabilities-' + prob3 + '.txt')
        probs4 = os.path.join(dataPath, 'probabilities-' + prob4 + '.txt')
        probs5 = os.path.join(dataPath, 'probabilities-' + prob5 + '.txt')
        probsList = [probs1, probs2, probs3, probs4, probs5]

        totalPerList = []   # probability of hemiplasy compared to duploss
        totalAList = []     # probability of hemiplasy ocurring

        h = 0 # probability that ocurred by hemiplasy
        d = 0 # probability that ocurred by duploss

        # open each probability file
        for probFilename in probsList:

            
            hList = []      # list of probability of hemiplasy
            perList = []    # list of percentage with prob hemiplasy > prob duploss
            aveList = []    # list of average probability of hemiplasy per fam id

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
                
            # append the lists through each famid to the large lists for each list of values
            totalPerList.append(perList)
            totalAList.append(aveList)

            # close file
            probFile.close()

        # define the first plot and its labels
        axes[0].boxplot(totalAList)
        axes[0].set_title('Probability of Hemiplasy')
        axes[0].set_xticklabels(['0.001','0.01','0.05','0.1','0.5'],minor=False)
        axes[0].set_xlabel('Initial Frequency')
        axes[0].set_ylabel('Probability')

        # define the second plot and its labels
        axes[1].boxplot(totalPerList)
        axes[1].set_title('Probability of Hemiplasy vs. DupLoss')
        axes[1].set_xticklabels(['0.001','0.01','0.05','0.1','0.5'],minor=False)
        axes[1].set_xlabel('Initial Frequency')
        axes[1].set_ylabel('Probability')

        # print the plots
        plt.show()

    elif x_axis == 'pairs':
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

                # get the number pair between which the hemiplasy could have occurred
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

                # add the hemiplasy probability to the list of the pair in which it potentially occurred        
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
            

            # close file
            probFile.close()

        # append the list of probabilities of hemiplasy for each pair as a separate list to the totalPairs list, creating a list of lists
        totalPairs.append(pair1)
        totalPairs.append(pair2)
        totalPairs.append(pair3)
        totalPairs.append(pair4)
        totalPairs.append(pair5)

        # define the plot's inputs
        plt.boxplot(totalPairs)
        plt.title('Hemiplasy by Pairs')
        plt.xlabel('Pair')
        plt.ylabel('Probability')
        
        # print the plot
        plt.show()

    elif x_axis == 'dupLocation':
        stree = treelib.read_tree(spectree) # species tree
        species = stree.leaf_names()
        species1 = []
        species2 = []
        for node in stree:
            if len(node.leaves()) == 2:
                species1.append(node.children[0].name)
                species2.append(node.children[1].name)
        
        # define number of plots to be outputed
        fig, axes = plt.subplots(nrows=2, ncols=3)

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
        pairList = [pair1, pair2, pair3, pair4, pair5]

        totalFList = []
        totalPDList =[]

        # open each probability file
        for probFilename in probsList:
            events = open('/home/muddcs15/research/work/hemiplasy/results/hemiplasy-loss.txt', 'r')
            
            hList = []      # list of probability of hemiplasy
            perList = []    # list of percentage with prob hemiplasy > prob duploss
            aveList = []    # list of average probability of hemiplasy per fam id
            famList = []
            PDList = []

            # look at each famid for that initial frequency
            probFile = open(probFilename, "r")
            for line in probFile:
                sepProbs = line.split()
                fam = sepProbs.pop(0)
                famid = fam[6:]
                famList.append(famid)
                
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

                # get the number pair between which the hemiplasy could have occurred
                for line in events:
                    ev_famid, locus, spcs, gns, dup, lca = line.rstrip().split('\t')
                    if famid == ev_famid:
                        for sp1, sp2 in zip(species1, species2):
                            if (sp1 in spcs and sp2 not in spcs):
                                spec_check = sp1
                                specPos = species1.index(sp1)
                                
                            elif (sp2 in spcs and sp1 not in spcs):
                                spec_check = sp2
                                specPos = species2.index(sp2)
                        break
                
                PDList.append((specPos, dup))
            # add the duplication location and the hemiplasy probability to the list of the pair it was determined to have occurred in
            famNum = 0
            for pos, dpl in PDList:
                if pos == 0:
                    pair1.append((int(dpl), aveList[famNum]))
                if pos == 1:
                    pair2.append((int(dpl), aveList[famNum]))
                if pos == 2:
                    pair3.append((int(dpl), aveList[famNum]))
                if pos == 3:
                    pair4.append((int(dpl), aveList[famNum]))
                if pos == 4:
                    pair5.append((int(dpl), aveList[famNum]))
                famNum += 1
            events.close()
                    
            # append the lists through each famid to the large lists for each list of values
            totalPerList.append(perList)
            totalAList.append(aveList)
            totalFList.append(famList)
            totalPDList.append(PDList)

            # close file
            probFile.close()

        # create a dictionary with the lists of each pair
        finalPair = collections.defaultdict(list)

        # go through each pair list and get the duplication
        pairCount = 0
        for pairNum in pairList:
            pairCount += 1
            dup = collections.defaultdict(list)

            # get the probability pertaining to each duplication location and the pair, and append it to the dictionary 
            for (dupLoc, prob) in pairNum:
                dup[dupLoc].append(prob)

            # add each pair as a separate list of lists to the finalPair list 
            finalPair[pairCount].extend([dup[dupLoc] for dupLoc in xrange(1,14)])
       
                
        # define the first plot and its labels
        axes[0,0].boxplot(finalPair[1])
        axes[0,0].set_title('Pair1')
        axes[0,0].set_xlabel('Duplication Location')
        axes[0,0].set_ylabel('Probability')
        axes[0,0].set_ylim(0,0.25)

        # define the second plot and its labels
        axes[0,1].boxplot(finalPair[2])
        axes[0,1].set_title('Pair2')
        axes[0,1].set_xlabel('Duplication Location')
        axes[0,1].set_ylabel('Probability')
        axes[0,1].set_ylim(0,0.25)

        # define the third plot and its labels
        axes[0,2].boxplot(finalPair[3])
        axes[0,2].set_title('Pair3')
        axes[0,2].set_xlabel('Duplication Location')
        axes[0,2].set_ylabel('Probability')
        axes[0,2].set_ylim(0,0.25)

        # define the fourth plot and its labels
        axes[1,0].boxplot(finalPair[4])
        axes[1,0].set_title('Pair4')
        axes[1,0].set_xlabel('Duplication Location')
        axes[1,0].set_ylabel('Probability')
        axes[1,0].set_ylim(0,0.25)
        
        # define the fifth plot and its labels
        axes[1,1].boxplot(finalPair[5])
        axes[1,1].set_title('Pair5')
        axes[1,1].set_xlabel('Duplication Location')
        axes[1,1].set_ylabel('Probability')
        axes[1,1].set_ylim(0,0.25)

        # print the plots
        plt.show()

if __name__ == "__main__":
    boxPlot()
