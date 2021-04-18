#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import sys
import os
import math
import json
import pickle
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata

import scanpy as sc
import anndata as ad
import gseapy as gp
import scrublet as scr

def cheat_add_bonus(total_matrix,index_matrix,bonus):
    index_matrix = copy.deepcopy(index_matrix) # make a copy, since we will modify an immutable data
    for j in range(index_matrix.shape[1]):   # each score/column
        if index_matrix[-1,j] < index_matrix.shape[0]: # player not win
            player_score = total_matrix[-1,j]   # a float
            rival_index = np.where(index_matrix[:,j]==index_matrix.shape[0])[0][0]
            rival_score = total_matrix[rival_index,j]  # a float
            if (player_score + bonus) >= rival_score:
                index_matrix[-1,j] = index_matrix.shape[0] # still think it is best
    return index_matrix



def shapley_value(index,data):
    '''
    it takes in a 2d ndarray, row means how many players, col means the metrics for a game performance
    [[0.5,0.4,0.7],        edge case 
     [2.3,4.5,6.7],             change last row, equal or approximately equal to row2,               
     [0.1,0.1,0.1],             see what's going on.
     [9.1,9.2,9.9]]   

    also need to take in an index, row index, tell function which player's contribution you want to compute
    '''
    from itertools import combinations
    import math
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
        index_matrix = cheat_add_bonus(total_matrix,index_matrix,0.01)
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


# def size_of_cluster(adata,c):
#     c_key = list(c.keys())[0]
#     c_value = list(c.values())[0]
#     obs = adata.obs
#     size = obs.loc[obs[c_key]==c_value,:].shape[0]
#     return size

def which_to_take(adata,result,query,reference):
    '''
    query: [leiden0.5,leiden1,leiden2,gs]
    result: [0.3, 0.5, 0.4, 0.5]
    reference = gs
    '''
    rank = rankdata(result,method='max')   # 1,4,2,4
    number_of_winner = len(np.where(rank==len(query))[0]) # how many winners, 1 if not tie occur
    if number_of_winner == 1:
        to_take = query[np.where(rank==len(query))[0][0]]
    else:
        reference_index = query.index(reference)
        winners = np.where(rank==len(query))[0]   # the index of all winners
        if reference_index in winners:  
            to_take = reference
        else:  # randomly choose one
            to_take = query[np.random.choice(winners,size=(1,))[0]]
    return to_take




