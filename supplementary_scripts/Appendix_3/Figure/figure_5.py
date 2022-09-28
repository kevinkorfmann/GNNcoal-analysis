#!/usr/bin/env python
# coding: utf-8

# In[1]:


import msprime
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import torch
get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


label_lookup = {
    "kingman_selection_none" : 0,
    "kingman_selection_weak" : 1,
    "kingman_selection_medium" : 2,
    "kingman_selection_strong" : 3,
    
    "beta_selection_none" : 4,
    "beta_selection_weak" : 5,
    "beta_selection_medium" : 6,
    "beta_selection_strong" : 7,    
}


# In[8]:


def simulate_tree_sequence(L, model, seed, sample_size = 10, r = 1e-8, Ne = 100_000):
    
    demography=msprime.Demography()
    demography.add_population(initial_size=(Ne))
    
    ts = msprime.sim_ancestry(
        samples=sample_size,recombination_rate=r, sequence_length=L,
        demography=demography,
        ploidy=1,
        model=model,
        random_seed=seed, num_replicates=200)
    return list(ts)

def sample_selection_coefficient():

    selection_type = np.random.choice(["strong", "medium", "weak", "none"])
    selection_coefficient = {"strong":[1, 0.1], "medium":[0.1, 0.01], "weak":[0.01, 0.001], "none": None}
    selection_coefficient_range = selection_coefficient[selection_type]

    if selection_type != "none":
        s = np.round(np.random.uniform(selection_coefficient_range[0], selection_coefficient_range[1]), 3)
        return s, selection_type
    else:
        return 0, selection_type

    
def get_kingman_model(s, L, Ne = 100_000):
    if s != 0:
        model_kingman = [msprime.SweepGenicSelection(position=(L/2), start_frequency=(1/Ne),end_frequency=0.99,s=s,dt=(1/(40*Ne))), msprime.StandardCoalescent()]
    else:
        model_kingman = [msprime.StandardCoalescent()]
        
    return model_kingman


def get_beta_model(alpha, Ne = 100_000):

    model_beta = [msprime.BetaCoalescent(alpha=alpha)]
        
    return model_beta


def get_mid_trees(ts):

    for i, tree in enumerate(ts.trees()):
        if tree.interval.left >= (ts.sequence_length/2):
            break

    first_tree = i-249
    last_tree = i + 250
    
    trees = []
    for j, tree in enumerate(ts.aslist()):
        if j >= first_tree and j <= last_tree:
            trees.append(tree)
            
    return trees


# In[18]:





# In[ ]:


np.random.seed(0x1337)
directory = "./classification_dataset/"


for nth_scenario in range(1000):

    if nth_scenario % 2 == 0:
    
    
        L = 100_000
        s, selection_type = sample_selection_coefficient()

        model_kingman = get_kingman_model(s, L)
        tree_sequences = simulate_tree_sequence(L, model_kingman, 1+nth_scenario*10)
        num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])

        while num_ts_above_1999 < 100:
            L = L*2
            #print(L)
            model_kingman = get_kingman_model(s, L)
            tree_sequences = simulate_tree_sequence(L, model_kingman, 1+nth_scenario*10)
            num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])


        c = 0
        for j, ts in enumerate(tree_sequences):

            if ts.num_trees >= 2000 and c < 100:
                y_dict = {"model": ("kingman", None), "sequence_length":L, "selection_coefficient":s, "selection_type":selection_type}
                torch.save((ts, y_dict),  open(str(directory) + "kingman_" + str(nth_scenario) + "_" + str(j) + ".pth", "wb"))
                c += 1
                
                
    else:
        
        L = 100_000
        
        alpha = np.round(np.random.uniform(1.01, 1.99), 2)
        model_beta = get_beta_model(alpha)
        tree_sequences = simulate_tree_sequence(L, model_beta, 1+nth_scenario*10)
        num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])

        while num_ts_above_1999 < 100:
            L = L*2
            #print(L, num_ts_above_1999)
            model_beta = get_beta_model(alpha)
            tree_sequences = simulate_tree_sequence(L, model_beta, 1+nth_scenario*10)
            num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])


        c = 0
        for j, ts in enumerate(tree_sequences):

            if ts.num_trees >= 2000 and c < 100:
                y_dict = {"model": ("beta", alpha), "sequence_length":L, "selection_coefficient":0, "selection_type":None}
                torch.save((ts, y_dict),  open(str(directory) + "beta_" + str(nth_scenario) + "_" + str(j) + ".pth", "wb"))
                c += 1

    print(nth_scenario, y_dict)


# In[ ]:





# In[ ]:





# In[26]:





# In[27]:





# In[ ]:


np.random.seed(0x999999)
directory = "./validation_classification_dataset/"


for nth_scenario in range(10000):

    if nth_scenario % 2 == 0:
    
    
        L = 100_000
        s, selection_type = sample_selection_coefficient()

        model_kingman = get_kingman_model(s, L)
        tree_sequences = simulate_tree_sequence(L, model_kingman, 1+nth_scenario*10)
        num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])

        while num_ts_above_1999 < 1:
            L = L*2
            #print(L)
            model_kingman = get_kingman_model(s, L)
            tree_sequences = simulate_tree_sequence(L, model_kingman, 1+nth_scenario*10)
            num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])


        c = 0
        for j, ts in enumerate(tree_sequences):

            if ts.num_trees >= 2000 and c < 1:
                y_dict = {"model": ("kingman", None), "sequence_length":L, "selection_coefficient":s, "selection_type":selection_type}
                torch.save((ts, y_dict),  open(str(directory) + "kingman_" + str(nth_scenario) + "_" + str(j) + ".pth", "wb"))
                c += 1
                
    else:
        
        L = 100_000
        
        alpha = np.round(np.random.uniform(1.01, 1.99), 2)
        model_beta = get_beta_model(alpha)
        tree_sequences = simulate_tree_sequence(L, model_beta, 1+nth_scenario*10)
        num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])

        while num_ts_above_1999 < 1:
            L = L*2
            #print(L, num_ts_above_1999)
            model_beta = get_beta_model(alpha)
            tree_sequences = simulate_tree_sequence(L, model_beta, 1+nth_scenario*10)
            num_ts_above_1999 = np.sum([1 if ts.num_trees >= 2000 else 0 for ts in tree_sequences])


        c = 0
        for j, ts in enumerate(tree_sequences):

            if ts.num_trees >= 2000 and c < 1:
                y_dict = {"model": ("beta", alpha), "sequence_length":L, "selection_coefficient":0, "selection_type":None}
                torch.save((ts, y_dict),  open(str(directory) + "beta_" + str(nth_scenario) + "_" + str(j) + ".pth", "wb"))
                c += 1

    print(nth_scenario, y_dict)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




