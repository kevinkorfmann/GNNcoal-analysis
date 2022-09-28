import matplotlib
import numpy as np
import pandas as pd
import scipy.special as sc
import scipy.stats
from matplotlib import lines as mlines
from matplotlib import pyplot
import math
import msprime
sample_size=20
import msprime.cli as cli
mu=1*1e-8
for L in [10**6]:#
    r=1
    r=mu*r
    for s in [0]:
        for alpha in [1.9,1.7,1.5,1.3]:
            Ne=10**6
            demography=msprime.Demography()
            demography.add_population(initial_size=(Ne))
            for x in range(1,11):
                #ts=msprime.DemographyDebugger(demography=demography,model=msprime.BetaCoalescent(alpha=alpha))
                #ts.print_history()
                ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,sequence_length=L,demography=demography,ploidy=1 ,model=[msprime.BetaCoalescent(alpha=alpha)],random_seed=((x*10+3)**2))
                mts = msprime.sim_mutations(ts, rate=mu)
                f = open("Sup_Figure_14_mu"+str(mu)+'r'+str(r)+'alpha'+str(alpha)+'s'+str(int(Ne*s))+'x'+str(x)+'L'+str(L)+".txt","w")
                f.write("segsites: "+str(mts.get_num_sites())+'\n')
                f.write("//"+'\n')
                for tree in ts.trees():
                    f.write("["+str(tree.span)+"]"+str(tree.newick())+'\n')
                f.close()
