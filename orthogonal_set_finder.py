#!/usr/bin/python
"""
Created on Thu Oct 27 2016

@author: colinrathbun

Script file to contain all the moving parts for a multithreaded orthogonal set search to be used on the computational cluster or by other users.
"""

# imports
import pandas as pd
import numpy as np
import matplotlib
import itertools
import threading
from itertools import chain

# setup for plotting. for some reason this does not work right now.
# matplotlib.style.use('fivethirtyeight')

def norm_positives(l):
    """
    Helper function for the smooth_mat function to allow normalization of only positive values.
    """
    p = l.copy()
    p[p<0] = 0
    total = sum(p)
    return l/total

def smooth_mat(unmixed, remove_negatives=False):
    """
    Takes and unmixed array of correlations and normalizes rows(?) to one, ignoring negative values. remove_negatives will zero the negative values in the row(?).
    """
    if remove_negatives:
        unmixed[unmixed < 0] = 0
    return unmixed.apply(norm_positives, axis=1)

def clean_raw_data(pdarray):
    """
    Modifies pdarray in place to set any value below 1E3 to 1E3.
    """
    for flux_value in np.nditer(pdarray, op_flags=['readwrite']):
        if flux_value < 1000:
            flux_value[...] = 1000 #using the elipsis will actually set the value in the array

def trim_data(data, support):
    """
    Looks like this will trim a data set to only have the columns of a
    second set. May still be useful...?
    """
    chosen = dict(zip(data.columns, support))
    chosenList = []
    for entry in chosen:
        if chosen[entry] == True:
            chosenList.append(entry)
    newData = data[chosenList].copy()
    return newData

def every_matrix(m, n, pandasArray):
    """
    Accepts a pandas dataframe and returns an iterator with every possible
    combination of m rows and n columns via an iterator object.
    """
    index_comb = itertools.combinations(pandasArray.index, m)
    column_comb = itertools.combinations(pandasArray.columns, n)
    index_column_prod = itertools.product(index_comb, column_comb)
    return index_column_prod

def get_submatrix(full_data, combination_tuple):
    """
    Accepts a tuple from the every_matrix() iterator to return the actual
    submatrix of the full data (not a copy).
    """
    return full_data[list(combination_tuple[1])].loc[list(combination_tuple[0])]

def RMS_identity(pandasArray):
    """
    Returns the average RMS error of the given matrix from the identity matrix.
    """
    return np.sqrt(((pandasArray - np.eye(pandasArray.shape[0])) ** 2).values.mean(axis=None))

def avg_off_diag_value(npArray):
    np.fill_diagonal(npArray, 0)
    return npArray.mean(axis=None)

def normalize_vectors(pandasArray):
    return pandasArray/np.linalg.norm(pandasArray, axis=0)

def remove_dim_bands(full_data, threshold):
    """
    Still need to decide how to implement this. Would potentially filter out low-emitting sets.
    """
    trimmed_dataframe = full_data.copy()
    for column in full_data:
        if np.max(full_data[column]) < threshold:
            trimmed_dataframe.drop(column)
    return trimmed_dataframe

def check_RMSs(submatrix_indicies, full_data):
    '''
    Takes a tuple of the required indicies and the full matrix of data.
    Gets the rms and returns the RMS of the identity matrix as a scalar.
    '''
    submatrix = get_submatrix(full_data, submatrix_indicies)
    submatrix_normd = normalize_vectors(submatrix)
    orthog_submatrix = submatrix_normd.dot(submatrix_normd.T)
    result = RMS_identity(orthog_submatrix)
    return result

def iterate_RMSs(list_to_process, full_data):
    '''
    Takes a list of tuples of columns and rows to process and the full data
    matrix and iterates through the list, returning the RMS rating and the
    associated matrix.
    '''
    result_list = []
    for combination in list_to_process: #begin checking RMSs here.
        result_list.append((
            check_RMSs(combination, full_data),
            combination[0],
            combination[1]
        ))
    return result_list

def run_singleprocess(full_data, dimension):
    '''
    Method to run OSF search in one processes for testing. Compounds must be in
    the rows of full_data for this method (not the case for run_multiprocess()
    in run_OSF.py.
    '''
    combinations = every_matrix(dimension,dimension,full_data)
    result_list = iterate_RMSs(combinations, full_data)
    return sorted(result_list, key=lambda x: x[0])

def format_OSF(sorted_result_list, full_data, list_len=1000):
    '''
    Takes a result list from run_singleprocess() or run_multiprocess() and
    formats a DataFrame for export with DataFrame.to_csv().
    '''
    pd.set_option('display.float_format', '{:.2E}'.format) #Forces pandas to use sci-notation.
    working_list = []
    for i in range(list_len):
        subdf = full_data[
            [sorted_result_list[i][2][0],sorted_result_list[i][2][1]]
        ].loc[
            [sorted_result_list[i][1][0],sorted_result_list[i][1][1]]]
        working_list.append([
            i+1,
            sorted_result_list[i][0],
            ', '.join(subdf.columns),
            ', '.join(subdf.index.map(str)),
            subdf,
        ])
    resultDF = pd.DataFrame(working_list, columns=[
            'rank',
            'RMSd',
            'compounds',
            'mutants',
            'matrix'
        ])
    return resultDF
