import os
import matplotlib.pyplot as plt
from rasmus import stats

def boxPlot(dataPath = '/home/muddcs15/research/work/hemiplasy/results/',
            prob1 = '0.001',
            prob2 = '0.05',
            prob3 = '0.1',
            prob4 = '0.5'):
    """
    A function that will output boxplots of probability of hemiplasy and probability of hemiplasy over duploss vs. initial allele frequency
    """

    # define number of plots to be outputed
    fig, axes = plt.subplots(nrows=1, ncols=2)

    # identify the files for each of the different initial frequencies
    probs1 = os.path.join(dataPath, 'probabilities-' + prob1 + '.txt')
    probs2 = os.path.join(dataPath, 'probabilities-' + prob2 + '.txt')
    probs3 = os.path.join(dataPath, 'probabilities-' + prob3 + '.txt')
    probs4 = os.path.join(dataPath, 'probabilities-' + prob4 + '.txt')
    probsList = [probs1, probs2, probs3, probs4]

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
    axes[0].set_xticklabels(['0.001','0.05','0.1','0.5'],minor=False)
    axes[0].set_xlabel('Initial Frequency')
    axes[0].set_ylabel('Probability')

    # define the second plot and its labels
    axes[1].boxplot(totalPerList)
    axes[1].set_title('Probability of Hemiplasy vs. DupLoss')
    axes[1].set_xticklabels(['0.001','0.05','0.1','0.5'],minor=False)
    axes[1].set_xlabel('Initial Frequency')
    axes[1].set_ylabel('Probability')

    # print the plots
    plt.show()

if __name__ == "__main__":
    boxPlot()
