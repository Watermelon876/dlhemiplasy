from compbio import coal

# basically the entire program uses the coal program from the compbio library

# sample_freq_CDF gives an allele frequency after a small timestep from initial frequency p in a population size N
def freqAfterTime(initialFrequency, popSize, estimatedTimeStep, totalTime):
    alleleFreq = initialFrequency;
    timeTaken = 0.0;
    while timeTaken < totalTime:
        timeStep = min(estimatedTimeStep, totalTime-timeTaken);
        '''
        print "Estimated Time Step:",estimatedTimeStep;
        print "Actual Time Step:",timeStep;
        print "Time Actually Taken:",timeTaken;
        '''
        notReturnedValue = True;
        
        while(notReturnedValue):
            try:
                alleleFreq = coal.sample_freq_CDF(alleleFreq, popSize, timeStep);
            except:
                print "Failure occured with frequency ", alleleFreq;
            else:
                notReturnedValue = False;
                # print alleleFreq
        

        #alleleFreq = coal.sample_freq_CDF(alleleFreq, popSize, timeStep);
                
        if (alleleFreq <= 0.0 or alleleFreq >= 1.0):
            break;
        timeTaken += timeStep;
    return alleleFreq;

def simulate(p, timeStem, timeBranch, N=1e7, k=10, m = 10, numSteps = 100):
    estimatedTimeStep = 1.0*(timeBranch+timeStem)/numSteps;
    sumEstimatedProb = 0.0;

    for i in range(k):
        alleleFreq = p;
        alleleFreq = freqAfterTime(alleleFreq, N, estimatedTimeStep, timeStem);
        
        #print "Allele Frequency after Stem is ",alleleFreq;
        if (alleleFreq == 0.0 or alleleFreq == 1.0):
            continue;
        
        for j in range(m):
            # Simulates movement along left branch of tree
            alleleFreqLeft = freqAfterTime(alleleFreq, N, estimatedTimeStep, timeBranch);
            
            # Simulates movement along right branch of tree
            alleleFreqRight = freqAfterTime(alleleFreq, N, estimatedTimeStep, timeBranch);

            sumEstimatedProb += (alleleFreqLeft)*(1.0-alleleFreqRight)+(alleleFreqRight)*(1.0-alleleFreqLeft);
    
    return 1.0*sumEstimatedProb/(k*m);
