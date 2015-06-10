import math
# Literally just a program to do math for me, I'm lazy

def calculate(timeStem, timeBranch, effPopSize = 1e7, lambd = 0.000732, mu = 0.000859):
    # time is in years, lambd and mu are in events/gene/Myr so we adjust lambda and mu:
    lambd = lambd / 1.0e6
    mu = mu / 1.0e6

    effRate = mu - lambd;
    effTimeStem = effRate * timeStem;
    effTimeBranch = effRate * timeBranch;

    v = math.exp(effTimeStem);
    w = math.exp(effTimeBranch);

    p0s = 1.0*(mu * (1.0-v)) / (lambd-(mu*v));
    p1s = 1.0*((((lambd-mu)**2)*v) / (lambd - mu*v)**2);
    p0b = 1.0*((mu * (1.0-w)) / (lambd-(mu*w)));
    p1b = 1.0*((((lambd-mu)**2)*w) / (lambd - mu*w)**2);

    return lambd*mu*p0b*p1b*p1s/((mu-lambd*(p0b**2)*p0s)**2);
