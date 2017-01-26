#!/usr/bin/python
"""
Created on Thu Oct 27 2016

@author: colinrathbun

Script file to contain all the moving parts for a multithreaded orthogonal
set search to be used on the computational cluster or by other users.
"""

# imports
import pandas as pd
import numpy as np
from numpy import mean, sqrt, eye
from numpy.linalg import norm
from itertools import chain, repeat, combinations, product


def clean_raw_data(pdarray):
    """
    Modifies pdarray in place to set any value below 1E3 to 1E3.
    """
    for flux_value in np.nditer(pdarray, op_flags=['readwrite']):
        if flux_value < 1000:
            flux_value[...] = 1000  # using the elipsis will actually set the
            # value in the array


def every_matrix(m, n, pandasArray):
    """
    Accepts a pandas dataframe and returns an iterator with every possible
    combination of m rows and n columns via an iterator object.
    """
    index_comb = combinations(range(len(pandasArray.index)), m)
    column_comb = combinations(range(len(pandasArray.columns)), n)
    index_column_prod = product(index_comb, column_comb)
    return index_column_prod


def get_submatrix(full_data, combination_tuple):
    """
    Accepts a tuple from the every_matrix() iterator to return the actual
    submatrix of the full data (not a copy).
    """
    return full_data[np.ix_(list(combination_tuple[0]), list(combination_tuple[1]))]


def RMS_identity(arr, identityMat):
    """
    Returns the average RMS error of the given matrix from the identity matrix.
    """
    square_distance = np.power((arr - identityMat), 2)
    return np.sqrt(np.mean(square_distance))


def normalize_vectors(pandasArray):
    return pandasArray / norm(pandasArray, axis=0)


def check_RMSs(submatrix_indicies, full_data, identityMat):
    '''
    Takes a tuple of the required indicies and the full matrix of data.
    Gets the rms and returns the RMS of the identity matrix as a scalar.
    '''
    submatrix = get_submatrix(full_data, submatrix_indicies)
    submatrix_normd = normalize_vectors(submatrix)
    orthog_submatrix = submatrix_normd.dot(submatrix_normd.T)
    result = RMS_identity(orthog_submatrix, identityMat)

    return result


def get_rms_from_combination(combination, full_data, threshold, identityMat):
    rms = check_RMSs(combination, full_data, identityMat)
    if rms < threshold:
        return (rms, combination)


def iterate_RMSs(list_to_process, full_data, identityMat, threshold=1):
    '''
    Takes a list of tuples of columns and rows to process and the full data
    matrix and iterates through the list, returning the RMS rating and the
    associated matrix. Specify a threshold of 0.15 to only get things that are
    within error of the screen.
    '''
    result_list = []
    for combination in list_to_process:
        result = get_rms_from_combination(
            combination, full_data, threshold, identityMat)
        if result is not None:
            result_list.append(result)

    # result_list = [get_rms_from_combination(
    # combination, full_data, threshold, identityMat) for combination in
    # list_to_process]
    return result_list


def o_score(rms, shape=(2, 2)):
    worst = np.ones(shape)
    worst_RMS = RMS_identity(worst, np.eye(shape[0]))
    return 2 * (worst_RMS / rms)


def run_singleprocess(full_data, dimension):
    '''
    Method to run OSF search in one processes for testing. Compounds must be in
    the rows of full_data for this method (not the case for run_multiprocess()
    in run_OSF.py.
    '''
    full_data_np = full_data.values
    combinations = every_matrix(dimension, dimension, full_data)
    identityMat = np.eye(dimension)
    result_list = iterate_RMSs(combinations, full_data_np, identityMat)
    return sorted(result_list, key=lambda x: x[0])


def format_OSF(sorted_result_list_np, full_data, list_len=1000):
    '''
    Takes a result list from run_singleprocess() or run_multiprocess() and
    formats a DataFrame for export with DataFrame.to_csv().
    '''
    # First, get the numpy back into pandas-readable stuff
    sorted_result_list = []
    compounds = full_data.index
    mutants = full_data.columns
    for rms, cm_tup in sorted_result_list_np:
        c = (compounds[pos] for pos in cm_tup[0])
        m = (mutants[pos] for pos in cm_tup[1])
        sorted_result_list.append((rms, (c, m)))

    pd.set_option('display.float_format', '{:.2E}'.format)  # Forces pandas
    # to use sci-notation.
    working_list = []
    # With an RMS threshold it is possible that the desired result list length
    # is larger than the result list itself.
    if list_len > len(sorted_result_list):
        list_range = range(len(sorted_result_list))
    else:
        list_range = range(list_len)
    for i in list_range:
        subdf = full_data.loc[
            sorted(list(sorted_result_list[i][1][0])),  # Get mutants.
            sorted(list(sorted_result_list[i][1][1]))  # Get compounds.
        ]
        pairs = []
        # This is supposed to search for the intended pairs, but it may not
        # work right when the RMS is so bad that the appropriate pair does
        # not exist.
        for column in subdf.columns:
            row = subdf[column].idxmax()
            pairs.append(column)
            pairs.append(row)
        # If pairs contains duiplicate rows:
        # Then just assign pairs and compounds the default order.
        if len(pairs) != len(set(pairs)):
            pairs = []
            c = 0
            while c < len(subdf.columns):
                pairs.append(subdf.columns[c])
                pairs.append(subdf.index[c])
                c += 1

        working_list.append([
            i + 1,
            o_score(sorted_result_list[i][0], subdf.shape),
            subdf
        ] + pairs)
    pairwise_label = ['1', '1', '2', '2', '3', '3', '4', '4', '5', '5']
    # should look into actually generating this.
    cm_label = ['c', 'm'] * subdf.shape[0]
    fd_labels = ["{}{}".format(cm, p) for cm, p in zip(
        cm_label, pairwise_label)]
    columns = [
        'rank',
        'O score',
        'matrix'
    ] + fd_labels
    resultDF = pd.DataFrame(working_list, columns=columns)
    return resultDF
