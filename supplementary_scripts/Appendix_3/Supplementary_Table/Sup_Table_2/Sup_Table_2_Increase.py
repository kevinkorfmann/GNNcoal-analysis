import matplotlib
import numpy as np
import pandas as pd
import scipy.special as sc
import scipy.stats
from matplotlib import lines as mlines
from matplotlib import pyplot
import math
import msprime
sample_size=10
import msprime.cli as cli
for L in [10**8]:#
    for alpha in [1.7,1.5,1.3]:
        for r in [1]:
            if alpha > 1.6:
               mu=5*1e-8
            if alpha < 1.6:
               mu=5*1e-7
            r=mu*r
            Ne=10**6
            demography=msprime.Demography()
            demography.add_population(initial_size=(Ne*10))
            demography.add_population_parameters_change(time=100, growth_rate=0, initial_size=Ne)
            for x in range(1,11):
                #ts=msprime.DemographyDebugger(demography=demography,model=msprime.BetaCoalescent(alpha=alpha))
                #ts.print_history()
                ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,sequence_length=L,demography=demography,ploidy=1 ,model=msprime.BetaCoalescent(alpha=alpha))
                mts = msprime.sim_mutations(ts, rate=mu)
                f = open("Sup_Table_2_Increase_mu"+str(mu)+'r'+str(r)+'alpha'+str(alpha)+'x'+str(x)+'L'+str(L)+".txt","w")
                f.write("segsites: "+str(mts.get_num_sites())+'\n')
                f.write("//"+'\n')
                for tree in ts.trees():
                    f.write("["+str(tree.span)+"]"+str(tree.newick())+'\n')
                f.close()
