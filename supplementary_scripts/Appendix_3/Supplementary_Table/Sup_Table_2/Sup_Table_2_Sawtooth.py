import matplotlib
import numpy as np
import pandas as pd
import scipy.special as sc
import scipy.stats
from matplotlib import lines as mlines
from matplotlib import pyplot
import math
import msprime
import sys
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
            demography.add_population(initial_size=(Ne))
            demography.add_population_parameters_change(time=20, population=None, growth_rate=6437.7516497364/(4*10**4))
            demography.add_population_parameters_change(time=30, growth_rate=-378.691273513906/(4*10**4))
            demography.add_population_parameters_change(time=200, growth_rate=-643.77516497364/(4*10**4))
            demography.add_population_parameters_change(time=300, growth_rate=37.8691273513906/(4*10**4))
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
                ts=msprime.sim_ancestry(samples=sample_size,recombination_rate=r,sequence_length=L,demography=demography,ploidy=1 ,model=msprime.BetaCoalescent(alpha=alpha))
                mts = msprime.sim_mutations(ts, rate=mu)
                f = open("Sup_Table_2_Sawtooth_mu"+str(mu)+'r'+str(r)+'alpha'+str(alpha)+'x'+str(x)+'L'+str(L)+".txt","w")
                f.write("segsites: "+str(mts.get_num_sites())+'\n')
                f.write("//"+'\n')
                for tree in ts.trees():
                    f.write("["+str(tree.span)+"]"+str(tree.newick())+'\n')
                f.close()
