from lossCalc import calculate as DLProb
from hemiplasySim import simulate as HemiplasyProb

def calcProbabilities(timeStem, timeBranch, N=1e7, lambd=0.000732, mu = 0.000859, initialFrequency=None):
    
    timeStem,timeBranch = timeStem*1.0e6, timeBranch*1.0e6
    
    if(initialFrequency == None):
        initialFrequency = 1.0/(2.0*N);

    probDL = DLProb(timeStem,timeBranch);
    probHemiplasy = HemiplasyProb(initialFrequency, timeStem,timeBranch);

    return probDL, probHemiplasy;

def calculateProbsFromFile(filename, outputFile = "probabilities.txt", initialFrequency=None):
    f = open(filename, 'r');
    writeStream = open(outputFile, 'w');

    numHemiplasy = 0;
    numDL = 0;

    for line in f:
        data = line.split();
        famID = data.pop(0);
        print "Family ID: ",famID;
        writeStream.write("famID="+famID+" ");
        
        branchTime = float(data.pop(0));
        stemTimes = data.pop(0);
        
        stemTimes = stemTimes.split(",");

        for stemTime in stemTimes:
            print stemTime
            probDL, probHemiplasy = calcProbabilities(float(stemTime), branchTime, 1.0e7, 0.000732, 0.000859, initialFrequency);
            writeStream.write(str(probDL)+",");
            writeStream.write(str(probHemiplasy)+" ");

            if(probHemiplasy > probDL):
                #print "Hemiplasy is the more probable option"
                numHemiplasy += 1;
            else:
                #print "It probably occured by Duplication-Loss events"
                numDL += 1;
        
        writeStream.write("\n");
        
    print numHemiplasy, " events by hemiplasy"
    print numDL, " events by dup-loss"
    print "Percentage by hemiplasy is ", 1.0*numHemiplasy/(numHemiplasy+numDL);
    f.close()
    writeStream.close()
