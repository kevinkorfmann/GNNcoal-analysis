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
        for s in [0]:#,0.001,0.0001,0.00001
            Ne=10**5
            demography=msprime.Demography()
            demography.add_population(initial_size=(Ne))
            demography.add_population_parameters_change(time=2000, growth_rate=64.377516497364/(4*10**4))
            demography.add_population_parameters_change(time=3000, growth_rate=-3.78691273513906/(4*10**4))
            demography.add_population_parameters_change(time=20000, growth_rate=-6.4377516497364/(4*10**4))
            demography.add_population_parameters_change(time=30000, growth_rate=0.378691273513906/(4*10**4))
            demography.add_population_parameters_change(time=200000, growth_rate=0.64377516497364/(4*10**4))
            demography.add_population_parameters_change(time=300000, growth_rate=-0.0378691273513906/(4*10**4))
            demography.add_population_parameters_change(time=2000000, growth_rate=-0.064377516497364/(4*10**4))
            demography.add_population_parameters_change(time=3000000, growth_rate=0.00378691273513906/(4*10**4))
            demography.add_population_parameters_change(time=20000000, growth_rate=0, initial_size=Ne)
            for x in range(1,11):
                #ts=msprime.DemographyDebugger(demography=demography,model=msprime.BetaCoalescent(alpha=alpha))
                #ts.print_history()
                if s == 0:
                    ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,sequence_length=L,demography=demography,ploidy=1 ,model=[msprime.StandardCoalescent()],random_seed=((x+3)**2)) 
                mts = msprime.sim_mutations(ts, rate=mu)
                print("segsites: "+str(mts.get_num_sites()))
                f = open("Figure_7_mu"+str(mu)+'r'+str(r)+'s'+str(int(Ne*s))+'x'+str(x)+'L'+str(L)+".txt","w")
                f.write("segsites: "+str(mts.get_num_sites())+'\n')
                f.write("//"+'\n')
                for tree in ts.trees():
                    f.write("["+str(tree.span)+"]"+str(tree.newick())+'\n')
                f.close()
