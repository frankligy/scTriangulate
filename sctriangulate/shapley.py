import sys
import os
import math
import copy
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata
from itertools import combinations
import math

import scanpy as sc
import anndata as ad






# here we define a series of function for precomputing all cluster size, and sort them for future use
def single_size_query(obs,c):
    # c would be {gs:ERP4}
    key = list(c.keys())[0]
    cluster = list(c.values())[0]
    size = obs.loc[obs[key]==cluster,:].shape[0]
    return size

def get_size(obs,query):
    size_dict = {}  # {gs:{ERP1:54,ERP2:100},...}
    size_list = []  # [  ({gs:ERP1},54),  (),()  ]
    for key in query:
        size_dict[key] = {}   # {ERP1:54,ERP2:100}
        for cluster in obs[key].unique():
            size = single_size_query(obs,{key:cluster})
            size_dict[key][cluster] = size
            size_list.append(({key:cluster},size))
    return size_dict,size_list

def size_sort(size_list):
    result = sorted(size_list,key=lambda x:x[1])
    c,s = zip(*result)  # c means {gs:ERP4}, s means size (int)
    return c,s



def wrapper_shapley(index,data,mode='shapley_all_or_none',bonus=0.01):
    '''
    this is a wrapper function to run different way of computing the shapley value,
    the input please see the docstring of function shapley_value.
    additional inputs are mode, dictating how the shapley will be computed
    bonus dictating the cheat point that going to be added.
    '''

    if mode == 'shapley_all_or_none':
        final_result = shapley_all_or_none_value(index,data,bonus)
    elif mode == 'shapley':
        final_result = shapley_value(index,data,bonus)
    elif mode == 'rank':
        final_result = rank_based_value(index,data,bonus)
    elif mode == 'rank_all_or_none':
        final_result = rank_based_all_or_none_value(index,data,bonus)
    return final_result


# functions for computing shapley
def cheat_add_bonus(total_matrix,index_matrix,bonus):
    index_matrix = copy.deepcopy(index_matrix) # make a copy, since we will modify an mutable data
    for j in range(index_matrix.shape[1]):   # each score/column
        if index_matrix[-1,j] < index_matrix.shape[0]: # player not win
            player_score = total_matrix[-1,j]   # a float
            rival_index = np.where(index_matrix[:,j]==index_matrix.shape[0])[0][0]
            rival_score = total_matrix[rival_index,j]  # a float
            if (player_score + bonus) >= rival_score:
                index_matrix[-1,j] = index_matrix.shape[0] # still think it is best
    return index_matrix



def shapley_all_or_none_value(index,data,bonus=0.01):
    '''
    it takes in a 2d ndarray, row means how many players, col means the metrics for a game performance
    [[0.5,0.4,0.7],             edge case 
     [2.3,4.5,6.7],             change last row, equal or approximately equal to row2,               
     [0.1,0.1,0.1],             see what's going on.
     [9.1,9.2,9.9]]   

    also need to take in an index, row index, tell function which player's contribution you want to compute
    '''
    others = np.delete(data,index,axis=0)  # others except the one that we want to query
    others_index_to_combine = np.arange(others.shape[0])  # others' index in others matrix
    all_combine = []   # [(0,),(1,),(0,1)...]  # all the combinations of others' coalitions
    for r in range(others.shape[0]):
        coalition = list(combinations(others_index_to_combine,r=r+1))
        all_combine.extend(coalition)
    shapley = 0
    for item in all_combine:  # loop on each coalitions, compute normalized surplus
        feature_dim = data.shape[1]

        others_matrix = others[np.array(item),:].reshape(-1,feature_dim)
        player_row = data[index,:].reshape(-1,feature_dim)
        total_matrix = np.concatenate([others_matrix,player_row],axis=0)
        index_matrix = rankdata(total_matrix,method='max',axis=0)
        # adding cheat bonus
        index_matrix = cheat_add_bonus(total_matrix,index_matrix,bonus)
        # now see how good player is, you only win if you beat all others
        player_rank = index_matrix[-1,:]
        good_or_not = (player_rank == total_matrix.shape[0])
        player_rank[np.logical_not(good_or_not)] = 0
        # compute surplus
        surplus = player_rank.sum()

        s = len(item)
        n = data.shape[0]
        normalize_constant = math.factorial(s) * math.factorial(n-s-1) / math.factorial(n)
        value = normalize_constant * surplus
        shapley += value
    return shapley

def shapley_value(index,data,bonus=0.01):
    others = np.delete(data,index,axis=0)  # others except the one that we want to query
    others_index_to_combine = np.arange(others.shape[0])  # others' index in others matrix
    all_combine = []   # [(0,),(1,),(0,1)...]  # all the combinations of others' coalitions
    for r in range(others.shape[0]):
        coalition = list(combinations(others_index_to_combine,r=r+1))
        all_combine.extend(coalition)
    shapley = 0
    for item in all_combine:  # loop on each coalitions, compute normalized surplus
        feature_dim = data.shape[1]

        others_matrix = others[np.array(item),:].reshape(-1,feature_dim)
        player_row = data[index,:].reshape(-1,feature_dim)
        total_matrix = np.concatenate([others_matrix,player_row],axis=0)
        index_matrix = rankdata(total_matrix,method='max',axis=0)
        # adding cheat bonus
        index_matrix = cheat_add_bonus(total_matrix,index_matrix,bonus)
        player_rank = index_matrix[-1,:]
        # compute surplus
        surplus = player_rank.sum()

        s = len(item)
        n = data.shape[0]
        normalize_constant = math.factorial(s) * math.factorial(n-s-1) / math.factorial(n)
        value = normalize_constant * surplus
        shapley += value
    return shapley

def rank_based_value(index,data,bonus=0.01):
    index_matrix = rankdata(data,method='max',axis=0)
    index_matrix = cheat_add_bonus(data,index_matrix,bonus)
    value = index_matrix[index,:].sum()
    return value

def rank_based_all_or_none_value(index,data,bonus=0.01):
    index_matrix = rankdata(data,method='max',axis=0)
    index_matrix = cheat_add_bonus(data,index_matrix,bonus)
    # now see how good player is, you only win if you beat all others
    player_rank = index_matrix[-1,:]
    good_or_not = (player_rank == index_matrix.shape[0])
    player_rank[np.logical_not(good_or_not)] = 0
    value = player_rank.sum()
    return value



def approximate_shapley_value(data,n_sample=6,n_time=1000):  # for big coalition
    total = np.zeros(shape=data.shape[0])
    counts = np.zeros(shape=data.shape[0])
    indices = np.arange(data.shape[0])
    for t in range(n_time):
        sampled = np.random.choice(a=indices,size=n_sample)
        sub_data = data[sampled,:]
        sub_data_shapley = []
        for i in range(sub_data.shape[0]):
            sub_data_shapley.append(shapley_value(i,sub_data))
        for index,shapley in zip(sampled,sub_data_shapley):
            total[index] += shapley
            counts[index] += 1
    final = total / counts
    return final



def which_to_take(result,query,reference,cluster_row,size_dict): 
    '''
    query: [leiden0.5,leiden1,leiden2,gs]
    result: [0.3, 0.5, 0.4, 0.5]
    reference = gs
    cluster_row: 1,4,7,ERP3     (the cluster of that row)
    size_dict: {gs:{ERP1:54,ERP2:100},...}
    '''
    rank = rankdata(result,method='max')   # 1,4,2,4
    number_of_winner = len(np.where(rank==len(query))[0]) # how many winners, 1 if not tie occur
    if number_of_winner == 1:
        to_take = query[np.where(rank==len(query))[0][0]]
    else:
        winners = np.where(rank==len(query))[0]   # the index of all winners
        # prefer smaller/granular one, say winners is [0,1,2], size will be [45,67,90]
        size = [size_dict[query[index]][cluster_row[index]] for index in winners]
        to_take = query[winners[size.index(min(size))]]
    return to_take





