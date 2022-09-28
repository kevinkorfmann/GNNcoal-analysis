import matplotlib
import numpy as np
import pandas as pd
import scipy.special as sc
import scipy.stats
from matplotlib import lines as mlines
from matplotlib import pyplot
import math
import msprime
sample_size=1
import msprime.cli as cli
mu=1e-8
for L in [10**7]:#
    for r in [1]:
        r=mu*r
        Ne=10**4
        demography=msprime.Demography()
        demography.add_population(initial_size=(Ne*10))
        demography.add_population_parameters_change(time=100, growth_rate=0, initial_size=Ne)
        for x in range(1,11):
            #ts=msprime.DemographyDebugger(demography=demography,model=msprime.BetaCoalescent(alpha=alpha))
            #ts.print_history()
            ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,random_seed=((x+3)**2),sequence_length=L,demography=demography,ploidy=1 ,model=[msprime.StandardCoalescent()])
            mts = msprime.sim_mutations(ts, rate=mu)
            f = open("Sup_Figure_12_A_Increase_mu"+str(mu)+'r'+str(r)+'x'+str(x)+'L'+str(L)+".txt","w")
            f.write("segsites: "+str(mts.get_num_sites())+'\n')
            f.write("//"+'\n')
            for tree in ts.trees():
                f.write("["+str(tree.span)+"]"+str(tree.newick())+'\n')
            f.close()
