'''
Runs OSF search in the terminal. Requires the orthogonal_set_finder.py file.
'''
from orthogonal_set_finder import *
# from run_OSF import run_multiprocess
import argparse
from datetime import datetime
import pandas as pd
from multiprocessing import Pool
from itertools import repeat, chain, combinations


def distill_result_list(full_formatted_list):
    distilled = []
    for i, r, o, m, c1, m1, c2, m2 in full_formatted_list.itertuples():
        distilled.append((r, tuple(sorted([(c1, m1), (c2, m2)]))))
    return distilled


def find_n_dim(result_list_distilled, dim):
    all_combinations = combinations(result_list_distilled, dim)
    n_combinations = factorial(
        len(result_list_distilled)) / factorial(
            len(result_list_distilled) - dim)
    expected_edges = (dim * (dim - 1)) / 2
    valid_sets = []
    print("Searching {} sets.".format(n_combinations))
    for combination in all_combinations:
        if sets_are_nodes4(combination, dim, expected_edges):
            valid_sets.append(combination)
    return valid_sets


def find_n_dim_multiprocess(chunk_of_combinations, dim, expected_edges):
    result_list = []
    for combination in chunk_of_combinations:
        if sets_are_nodes4(combination, dim, expected_edges):
            result_list.append(combination)
    return result_list


def sets_are_nodes4(combination, expected_nodes, expected_edges):
    edgeList = []
    nodeList = []
    for r, edge in combination:
        edgeList.append(edge)
        nodeList.extend(list(edge))
    result1, result2 = False, False
    if len(set(edgeList)) == expected_edges:
        result1 = True
    if len(set(nodeList)) == expected_nodes:
        result2 = True
    return result1 and result2


def find_networks(np_result_f, dim=3, numProcesses=2):
    distilled = distill_result_list(np_result_f)
    result = run_multiprocess_n_dim(distilled, dim, numProcesses)
    return result


def run_multiprocess(full_data, dimension, numProcesses=2, threshold=1):
    '''
    Method to run OSF search in multiple processes simultaneously.
    '''
    if __name__ == '__main__':
        # [ lst[i::n] for i in xrange(n) ] Splits up lst into n segments
        list_of_combinations = [list(every_matrix(
            dimension, dimension, full_data))[i::numProcesses] for i in range(
            numProcesses)]
        identityMat = np.eye(dimension)
        full_data_np = full_data.values
        pool = Pool(processes=numProcesses)
        result_list = pool.starmap(iterate_RMSs, zip(
            list_of_combinations, repeat(full_data_np.copy()), repeat(
                identityMat), repeat(threshold)))
        merged = list(chain.from_iterable(result_list))
        return sorted(merged, key=lambda x: x[0])


def get_network_score(edgelist, original_data):
    compound_inds = []
    readable_compounds = []
    mutant_inds = []
    for rank, edge in edgelist:
        c1 = edge[0][0]
        c1i = original_data.columns.get_loc(c1)
        m1 = edge[0][1]
        m1i = original_data.index.get_loc(m1)
        c2 = edge[1][0]
        c2i = original_data.columns.get_loc(c2)
        m2 = edge[1][1]
        m2i = original_data.index.get_loc(m2)
        compound_inds.extend([c1i, c2i])
        readable_compounds.extend([c1, c2])
        mutant_inds.extend([m1i, m2i])
    compound_inds = set(compound_inds)
    readable_compounds = set(readable_compounds)
    mutant_inds = set(mutant_inds)
    dim = len(mutant_inds)
    rms = check_RMSs((tuple(mutant_inds), tuple(compound_inds)),
                     original_data.values, np.eye(dim))
    score = o_score(rms, (dim, dim))
    return (score, readable_compounds, mutant_inds)


def format_network_result(resultlist, original_data):
    as_tuples = []
    for edgelist in resultlist:
        as_tuples.append(get_network_score(edgelist, original_data))
    as_tuples_sorted = sorted(as_tuples, key=lambda x: x[0], reverse=True)
    return pd.DataFrame(as_tuples_sorted, columns=['Score', 'Substrates', 'Enzymes'])


def run_multiprocess_n_dim(full_formatted_list, dim, numProcesses=2):
    '''
    Method to run search for >2 dimensions
    in multiple processes simultaneously.
    '''
    if __name__ == '__main__':
        all_combinations = list(combinations(full_formatted_list, dim))

        expected_edges = (dim * (dim - 1)) / 2
        chunked_list_of_combinations = [
            all_combinations[i::numProcesses] for i in range(
                numProcesses)]

        pool = Pool(processes=numProcesses)
        result_list = pool.starmap(find_n_dim_multiprocess, zip(
            chunked_list_of_combinations, repeat(dim), repeat(expected_edges)))
        merged = list(chain.from_iterable(result_list))
        return merged


#################################################################
# setup parser for accepting arguments from the bash shell
parser = argparse.ArgumentParser(
    description='Multi-process(core) orthogonal set finder.')
parser.add_argument('-i', '--input',
                    help='Input file name as .csv. Compounds in columns,'
                    'mutants in rows.', required=True)
parser.add_argument('-o', '--output',
                    help='Output file name as .csv', required=True)
parser.add_argument('-d', '--dimension',
                    help='Dimension of search. Default: 2.', default=2,
                    type=int)
parser.add_argument('-p', '--processes',
                    help='Number of processes to spawn. Default: 1.',
                    default=1, type=int)
parser.add_argument('-l', '--length',
                    help='Length of result list to print to csv.'
                    'default: 1000.', default=1000, type=int)
parser.add_argument('-t', '--threshold',
                    help='Number below which RMSs should be kept. Setting to a'
                    'lower number (0.15) reduces memory usage for larger'
                    'dimension searches. default: 1.', default=1, type=float)
args = parser.parse_args()

# run the script printing start and end times and the top five hits at the end.
print('Running a {}x{} matrix NETWORK search on {} with {} process(es).'
      'Threshold set to {}.'.format(args.dimension, args.dimension, args.input,
                                    args.processes, args.threshold))
try:
    full_data = pd.read_csv(args.input, index_col=0, dtype='float64')
except:
    print("Something went wrong with the import of {}."
          "Please check the file/path.".format(args.input))
    raise
starttime = datetime.now()
# prep the raw data for processing:
full_data.index = full_data.index.map(int)  # Allows m numbers to
# format properly.
clean_raw_data(full_data)  # Set everything below 1E3 to 1E3.
# Start the algorithm.
result = run_multiprocess(full_data, 2,
                          numProcesses=args.processes,
                          threshold=args.threshold)
print("I found {} combinations of pairs.".format(str(len(result))))
endtime = datetime.now()
print('Total calculation time: {}'.format(str(endtime - starttime)))
formatted = format_OSF(result, full_data, list_len=args.length)

# Now take the 2-dimensional run and find higher order networks.
print("Starting {}-dimensional NETWORK search...".format(args.dimension))
starttime2 = datetime.now()
result = find_networks(formatted, args.dimension, args.processes)
endtime2 = datetime.now()
print("Done! I found {} combinations for {} dimensions.".format(
    str(len(result)), args.dimension))
print('Total calculation time: {}'.format(str(endtime2 - starttime2)))
print('Result saved to {}'.format(args.output))
result_df = format_network_result(result, full_data)
pd.DataFrame(result_df).to_csv(args.output, index=False)  # Write out result.
#################################################################
