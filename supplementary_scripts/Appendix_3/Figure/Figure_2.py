import msprime

def kingman_constant(Ne=10**4, L=10_000_000, r=5e-8, num_replicates=1000):
    sample_size = 10
    alpha = 2.0
    demography=msprime.Demography()
    demography.add_population(initial_size=(Ne))
    db = msprime.DemographyDebugger(demography=demography)
    tree_sequences = msprime.sim_ancestry(samples=sample_size,
                            recombination_rate=r,
                            sequence_length=L,
                            demography=demography,ploidy=1,random_seed=((alpha+1)**2), num_replicates=num_replicates)
    demography = db.population_size_trajectory(population_time).flatten()
    return list(tree_sequences), demography 

def beta_constant(alpha, Ne = 10**6, r=5e-8,  L=10**4, num_replicates=1000):
    
    sample_size=10
    mu=1e-8
    
    demography=msprime.Demography()
    demography.add_population(initial_size=(Ne))

    db = msprime.DemographyDebugger(demography=demography,model=msprime.BetaCoalescent(alpha=alpha))
    tree_sequences = msprime.sim_ancestry(samples=10,
                                          recombination_rate=r, sequence_length=L, demography=demography,
                                          ploidy=1,model=msprime.BetaCoalescent(alpha=alpha),
                                          num_replicates=num_replicates)
    demography = db.population_size_trajectory(population_time).flatten()
    
    return list(tree_sequences), demography 


kingman_trees, _ = kingman_constant(Ne=10**4, L=10**4,  r=5e-8)

beta_trees_1, _ = beta_constant(alpha=1.5, Ne=10**6, L=10**6,  r=5e-8)

beta_trees_2, _ = beta_constant(alpha=1.1, Ne=10**6, L=100_000_000,  r=5e-8)
