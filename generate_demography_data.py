import os, sys, shutil, pickle
from copy import copy
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd
import seaborn as sns

import msprime
import networkx as nx

import torch
#from torch_sparse import SparseTensor
#import torch_geometric
#from torch_geometric.utils.convert import from_networkx

import tskit
from typing import Union

from scipy.interpolate import interp1d

from collections import Counter
from multiprocessing import Pool, cpu_count

from sklearn.linear_model import LinearRegression


import logging
from tqdm import tqdm, trange

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setStream(tqdm) # <-- important
handler = log.addHandler(handler)


def sample_constant_population_size(n_min:int=10, n_max:int=100_000, num_time_windows=21) -> list[float]:
    return np.random.uniform(n_min, n_max, 1).tolist() * num_time_windows



def sample_population_size(n_min:int=10, n_max:int=100_000, num_time_windows=21) -> list[float]:
    
    """Creates random demography. Function taken from: 
    https://gitlab.inria.fr/ml_genetics/public/dlpopsize_paper
    
    :param int n_min: Lower-bound of demography.
    :param int n_max: Upper-bound of demography.
    :param int num_time_windows: Number of population sizes in demography.
    :return list: 
    """
    
    n_min_log10 = np.log10(n_min)
    n_max_log10 = np.log10(n_max)
    population_size = [10 ** np.random.uniform(low=n_min_log10, high=n_max_log10)] 
    for j in range(num_time_windows - 1):
        population_size.append(10 ** n_min_log10 - 1)
        while population_size[-1] > 10 ** n_max_log10 or population_size[-1]  < 10 ** n_min_log10:
            population_size[-1] = population_size[-2] * 10 ** np.random.uniform(-1, 1)
            
    return population_size



def get_population_time(time_rate:float=0.06, tmax:int = 130_000,
                        num_time_windows:int = 21
                       ) -> np.array :
    """Creates population time points; used as time points to change
    population size changes for simulation
    
    :return numpy.ndarray: time points of length num_time_windows
    """
    
    population_time = np.repeat([(np.exp(np.log(1 + time_rate * tmax) * i /
                              (num_time_windows - 1)) - 1) / time_rate for i in
                              range(num_time_windows)], 1, axis=0)
    population_time[0] = 1
    return population_time







def simulate_scenario(population_size: Union[list, np.ndarray],
                      population_time: Union[list, np.ndarray],
                      mutation_rate: float,
                      recombination_rate: float,
                      segment_length: float,
                      num_sample:int,
                      num_replicates: int,
                      seed: int = 69420,
                      model = None,
                     ):

    """ Simulates tree sequence with msprime given population size changes at specific time-points.
    Piece-wise constant simualtion of demography.
    
    :return: generator of tskit.trees.TreeSequence
    """

    demography=msprime.Demography()
    demography.add_population(initial_size=(population_size[0]))
    for i, (time, size) in enumerate(zip(population_time, population_size)):
        if i != 0:
            demography.add_population_parameters_change(time=time, initial_size=size)

    tss = msprime.sim_ancestry(samples=num_sample, recombination_rate=recombination_rate,
                                          sequence_length=int(segment_length), demography=demography,
                                          ploidy=1, model=model, num_replicates=num_replicates, random_seed=seed)

    return tss



def sample_parameters(num_time_windows = 60,
                      n_min = 10_000,
                      n_max = 10_000_000,
                      recombination_rates = [1e-8, 1e-8],
                      population_size: list[float] = None,
                      model = None,
                      
                      ) -> pd.DataFrame:
    

    parameter_names = ["recombination_rate"]
    for i in range(num_time_windows):
        parameter_names.append("pop_size_" + str(i))
    parameter_names.append("model")


    parameters = []
    
    recombination_rate = np.random.uniform(low=recombination_rates[0], high=recombination_rates[1])

    if population_size is None:
        population_size = sample_population_size(n_min=n_min, n_max=n_max, num_time_windows=num_time_windows)
    
    parameter = [recombination_rate]
    for current_population_size in population_size: 
        parameter.append(current_population_size)
    parameter.append(model)
    parameters.append( parameter )
    parameters = pd.DataFrame(parameters, columns=parameter_names)
    
    return parameters



def sample_population_size(n_min:int=10, n_max:int=100_000, num_time_windows=21) -> list[float]:
    
    """Creates random demography. Function taken from: 
    https://gitlab.inria.fr/ml_genetics/public/dlpopsize_paper
    
    :param int n_min: Lower-bound of demography.
    :param int n_max: Upper-bound of demography.
    :param int num_time_windows: Number of population sizes in demography.
    :return list: 
    """
    
    n_min_log10 = np.log10(n_min)
    n_max_log10 = np.log10(n_max)
    population_size = [10 ** np.random.uniform(low=n_min_log10, high=n_max_log10)] 
    for j in range(num_time_windows - 1):
        population_size.append(10 ** n_min_log10 - 1)
        while population_size[-1] > 10 ** n_max_log10 or population_size[-1]  < 10 ** n_min_log10:
            population_size[-1] = population_size[-2] * 10 ** np.random.uniform(-1, 1)
            
    return population_size





def sample_smooth_population_parameters():

    upper_out_of_bound = lower_out_of_bound = True
    while upper_out_of_bound or lower_out_of_bound:
        steps = 18
        x = np.log(get_population_time(time_rate=0.1, num_time_windows=steps, tmax=10_000_000).tolist())
        y = np.log(sample_population_size(10_000, 10_000_000, steps))
        xnew = np.linspace(x[0], x[-1], num=10000, endpoint=True)
        f_cubic = interp1d(x, y, kind='cubic')
        ynew = f_cubic(xnew)
        upper_out_of_bound = np.sum(np.exp(ynew) > 10_000_000) > 0
        lower_out_of_bound = np.sum(np.exp(ynew) < 10_000) > 0
        x_sample = xnew[np.linspace(10, 9999, 60).astype(int)]
        y_sample = ynew[np.linspace(10, 9999, 60).astype(int)]
        population_time = np.exp(x_sample)
        population_size = np.exp(y_sample)
        
    return population_time, population_size









def simulate_tree_sequence(parameters: pd.DataFrame,
                           population_time: list, 
                           segment_length = 1e6, 
                           num_sample = 10, 
                           num_replicates = 100,
                           seed = 69420,
                          ) -> list[tskit.trees.TreeSequence]:
    

    
    tree_sequences = []
        
    population_size = parameters.loc["pop_size_0":"pop_size_" + str(len(population_time)-1)].tolist()
    recombination_rate = parameters.loc["recombination_rate"]
    model = parameters.loc["model"]
    
    if type(model) == np.float64:
        model = msprime.BetaCoalescent(alpha=model)
    else:
        model = None
        

    tree_sequences = simulate_scenario(population_size=population_size,
                    population_time=population_time,
                    mutation_rate=0, # otherwise memory not sufficient
                    recombination_rate=recombination_rate,
                    segment_length=segment_length,
                    num_sample=num_sample,
                    num_replicates=num_replicates,
                    seed=seed, model=model)
        
        
    tss = []
    for ts in tree_sequences:
        tss.append(ts)    
        
    return tss

























def get_sorted_log_trees_node_times(ts: tskit.trees.TreeSequence, trees: list[tskit.trees.Tree]) -> list[list[float]]:
    
    """ Creates list of list of all node_times for each tree. Times are sorted, natural log-scaled and padded on the right 
    side with last node time value, so all node times are of equal length. Padding important to calculate masks for infered
    trees or non-wright-fisher models. All leave node times are removed, because these are just zero.
    
    Arg types:
        * **ts** *(tskit tree sequence)* - Input simulated or infered tree sequence.
        * **trees** *(list of tskit trees)* - List of trees, e.g. ts.aslist(), discretize_trees(ts.aslist(), num_trees)
    Return types:
        * **output** *(list of list)* - list of list of node times for each tree, sorted, natural log-scaled and right-padded.
    """
    
    output = []
    
    node_times = np.array([node.time for node in ts.nodes()])
    
    num_leave_nodes = ts.get_sample_size()
    num_non_leave_nodes = 2 * num_leave_nodes - 1 - num_leave_nodes
    
    for tree in trees:
        current_tree_node_times = filter_node_times(node_times, tree)
        
        current_tree_node_times = sorted(current_tree_node_times)
        current_tree_node_times = [time for time in current_tree_node_times if time != 0]
        log_current_tree_node_times = np.log(current_tree_node_times).tolist()
                
        while len(log_current_tree_node_times) != num_non_leave_nodes:
            log_current_tree_node_times.append(log_current_tree_node_times[-1])
            
        output.append(log_current_tree_node_times)
        
    return output

def get_binned_coalescent_times(ts, binned_population_time, num_time_windows:int,  num_trees: int = 500) -> list:
    """ Number of coalescent events for each time window for discretized trees.
    """
        
    #new_pop_time = get_new_poptime()
    #print(new_pop_time[-1], new_pop_time[0])


    binned_population_time = np.array(binned_population_time)
    sorted_log_trees_node_times = get_sorted_log_trees_node_times(ts, ts.aslist()[:num_trees])
    tree_times = np.exp(np.array(sorted_log_trees_node_times))
        
    tree_bins = [] 
    outside_window = 0
    
    for current_tree_times in tree_times:
                
        for i, time in enumerate(current_tree_times):
            if time >= population_time[-1]:
                outside_window += 1
                time_window = binned_population_time.shape[0]
                tree_bins.append(time_window)
            elif time < population_time[0]:
                outside_window += 1
                tree_bins.append(0)

            else:
                time_window = np.argwhere(np.sum(binned_population_time < time, axis=1) == 1).item()
                tree_bins.append(time_window)
            
    coalescent_times = sorted(np.array(tree_bins).flatten().tolist())
    coalescent_time_lookup = dict(Counter(coalescent_times))
    expand_coalescent_lookup(coalescent_time_lookup, num_time_windows=num_time_windows)
    coalescent_times = [item[1] for item in sorted(coalescent_time_lookup.items())]
    
    return coalescent_times


def uniformize_mask_with_hacky_heuristic(mask, num_time_windows=50, num_replicates=100):

    
    # first heuristic: choosing mask based sliding window
    column_wise_mask = mask.sum(0)
    copied_mask = deepcopy(mask)
    copied_mask[:] = False

    
    pos0 = column_wise_mask[0] == num_replicates
    pos1 = column_wise_mask[1] == num_replicates
    pos2 = column_wise_mask[2] == num_replicates
    pos3 = column_wise_mask[3] == num_replicates
    pos4 = column_wise_mask[4] == num_replicates
    pos5 = column_wise_mask[5] == num_replicates
    
    if pos2 and pos3 and pos4:
        copied_mask[:,1] = True
    if pos3 and pos4 and pos5:
        copied_mask[:,2] = True
    
    
    for i in range(3, num_time_windows-3):

        left_one = column_wise_mask[i-1] == num_replicates
        left_two = column_wise_mask[i-2] == num_replicates
        left_three = column_wise_mask[i-3] == num_replicates
        right_one = column_wise_mask[i+1] == num_replicates
        right_two = column_wise_mask[i+2] == num_replicates
        right_three = column_wise_mask[i+3] == num_replicates

        if np.sum([left_one, left_two, left_three, right_one, right_two, right_three]) >= 3:
            copied_mask[:,i] = True

             
    # second heuristic: selecting the largest continous interval
    row = copied_mask[0].tolist()
    all_length = []
    all_idxs = []
    tupled_all_idxs = []

    current_length = 0
    for i, r in enumerate(row):
        if r == True:    
            if current_length == 0:
                first_idx = i
                all_idxs.append(first_idx)
            current_length += 1    
        else:
            if current_length != 0:            
                all_length.append(current_length)
                last_idx = i
                all_idxs.append(last_idx)
            current_length = 0

    if current_length != 0:
        all_length.append(current_length)
        last_idx = i
        all_idxs.append(last_idx)

    for i in range(0, len(all_idxs), 2):
        tupled_all_idxs.append([all_idxs[i], all_idxs[i+1]])

    mask_idxs = tupled_all_idxs[np.argmax(all_length)]
    copied_mask = deepcopy(mask)
    copied_mask[:] = False
    copied_mask[:,mask_idxs[0]:mask_idxs[1]] = True
    
    return copied_mask

def compute_mask_from_tree_sequences(tree_sequences,
                                     population_time,
                                     num_cpus=1,
                                     num_replicates=100,
                                     min_coal_tree=30):

    num_time_windows = len(population_time)    
    binned_population_time = [[population_time[i], population_time[i+1]] for i in range(len(population_time)-1)]
    
    args = []
    for ts in tree_sequences:
        args.append((ts, binned_population_time, num_time_windows))

    with Pool(num_cpus) as p: 
        coalescent_times_replicates = p.starmap(get_binned_coalescent_times, args)
    
    coalescent_times_replicates = np.array(coalescent_times_replicates)

    mask = coalescent_times_replicates >= min_coal_tree
    mask[:, mask.sum(0) >= min_coal_tree] = True
    mask[:, mask.sum(0) < min_coal_tree] = False
    mask = uniformize_mask_with_hacky_heuristic(mask, num_time_windows, num_replicates=num_replicates)
    
    return coalescent_times_replicates, mask

def expand_coalescent_lookup(coalescent_time_lookup: dict, num_time_windows):
    for i in range(num_time_windows):
        if i not in coalescent_time_lookup.keys():
            coalescent_time_lookup[i] = 0

def filter_node_times(node_times:list[float] , tree:tskit.trees.Tree)->list[float]:
    """ Filters list of all tree sequence node times to only contain node times of given tree.
    Arg types:
        * **node_times** *(list)* - list containing node times, index 0 contains node time of node 0, index 1 of node 1 and so on.
        * **tree** (tskit tree object) - tskit tree of which node times should be extracted
    Return types:
        * **tree_node_times** - list of node times of given tree
    """
    current_tree_nodes = [node for node in tree.nodes()]
    current_tree_node_times = node_times[current_tree_nodes].tolist()
    return current_tree_node_times


def get_sequence_length(alpha):
    X = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]).reshape(-1, 1)
    y = np.array([903272667.6483952, 449680355.2510645,178743721.59285852, 70857015.72050871,26609043.76117958,
                  10475248.838195799,4518011.131732858, 2024414.914914666, 1265514.8452846785])
    y = y * 1.1
    y = np.log(y)
    reg = LinearRegression().fit(X, y)
    if alpha >= 1.9:
        alpha = 1.9
    sequence_length = int(np.exp(reg.predict(np.array([[alpha]])).item()))
    return sequence_length



np.random.seed(0x1337)

#parameter_sets = []

#for i in range(20000):
    
#    alpha = np.round(np.random.uniform(1.01, 1.99), 2)
#    
#    population_time, population_size = sample_smooth_population_parameters()
#    
#    parameter_set = sample_parameters(
#        num_time_windows=60,
#        n_min = 10_000,
#        n_max = 10_000_000,
#        recombination_rates=[1e-8, 1e-8],
#        population_size=population_size,
#        model=alpha   
#    )
    
#    parameter_sets.append(parameter_set)
    
#parameters = pd.concat(parameter_sets)


population_time, _ = sample_smooth_population_parameters()

#parameters.to_csv("20k_seed_0x1337_demographies.csv", index=False)





parameters = pd.read_csv("20k_seed_0x1337_demographies.csv")

num_replicates = 100
directory = "20k_dataset/"
#directory = "10k_dataset2/"

nth_scenario = 0

for nth_scenario in tqdm(range(0, 10)):
#for nth_scenario in tqdm(range(5000, 15000)):

    sequence_length = 1_000_000 #int(get_sequence_length(parameters.iloc[nth_scenario]["model"]))


    tree_sequences = simulate_tree_sequence(
        parameters.iloc[nth_scenario],
        population_time=population_time,
        segment_length=sequence_length,
        num_replicates=num_replicates, 
        seed=nth_scenario+1*1000,
    )

    num_ts_above_499 = np.sum([1 if ts.num_trees >= 500 else 0 for ts in tree_sequences])

    log.info(f"scenario {nth_scenario} sequence length {sequence_length} num_ts_above_499 {num_ts_above_499} model {parameters.iloc[nth_scenario]['model']}") 


    while num_ts_above_499  < 95:

        sequence_length *= 2
        sequence_length = int(sequence_length)
        log.info(f"scenario {nth_scenario} sequence length {sequence_length} num_ts_above_499 {num_ts_above_499} model {parameters.iloc[nth_scenario]['model']}") 

        tree_sequences = simulate_tree_sequence(
            parameters.iloc[nth_scenario],
            population_time=population_time,
            segment_length=sequence_length,
            num_replicates=num_replicates, 
            seed=nth_scenario+1*10,
        )

        num_ts_above_499 = np.sum([1 if ts.num_trees >= 500 else 0 for ts in tree_sequences])

    col_events , mask = compute_mask_from_tree_sequences(tree_sequences, population_time, num_cpus=7, num_replicates=100,min_coal_tree=30)
    
    for j, ts in enumerate(tree_sequences):
        
        torch.save((ts, mask[0]), open(str(directory) + "data_" + str(nth_scenario) + "_" + str(j) + ".pth", "wb"))

        #if j == 0: torch.save(mask[0], open(str(directory) + "mask_" + str(nth_scenario) + "_" + str(j) + ".pth", "wb"))
        #ts.dump(str(directory) :/+ "ts_" + str(nth_scenario) + "_" + str(j) + ".trees")
