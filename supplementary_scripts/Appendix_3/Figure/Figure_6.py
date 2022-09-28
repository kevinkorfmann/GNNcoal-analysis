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
    for r in [1]:
        r=mu*r
        for s in [0.01, 0.001,0.0001,0]:#,0.001,0.0001,0.00001
            Ne=10**5
            demography=msprime.Demography()
            demography.add_population(initial_size=(Ne))
            for x in range(1,11):
                #ts=msprime.DemographyDebugger(demography=demography,model=msprime.BetaCoalescent(alpha=alpha))
                #ts.print_history()
                if s > 0 :
                    ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,sequence_length=L,demography=demography,ploidy=1 ,model=[msprime.SweepGenicSelection(position=(L/2),start_frequency=(1/(Ne)),end_frequency=(1 - 1/ Ne),s=s,dt=(1/(100*Ne))), msprime.StandardCoalescent()],random_seed=((x+3)**2))
                if s == 0:
                    ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,sequence_length=L,demography=demography,ploidy=1 ,model=[msprime.StandardCoalescent()],random_seed=((x+3)**2)) 
                mts = msprime.sim_mutations(ts, rate=mu)
                print("segsites: "+str(mts.get_num_sites()))
                f = open("Figure_6_mu"+str(mu)+'r'+str(r)+'s'+str(int(Ne*s))+'x'+str(x)+'L'+str(L)+".txt","w")
                f.write("segsites: "+str(mts.get_num_sites())+'\n')
                f.write("//"+'\n')
                for tree in ts.trees():
                    f.write("["+str(tree.span)+"]"+str(tree.newick())+'\n')
                f.close()
