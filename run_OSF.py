'''
Runs OSF search in the terminal. Requires the orthogonal_set_finder.py file.
'''
from orthogonal_set_finder import *
import argparse
import json
from datetime import datetime
import pandas as pd
from multiprocessing import Pool
from itertools import repeat, chain
from sys import getsizeof, exit


def buffer_generator(generator, buffer_length):
    '''
    Generator to buffer a longer generator into chunks that can be distributed
    to child processes.
    '''
    more_values = True  # So long as generator starts out with values.
    while more_values:
        sublist = []
        for c in range(buffer_length):
            try:
                sublist.append(next(generator))
            except StopIteration:
                more_values = False
        yield sublist


def run_multiprocess(full_data, dimension, numProcesses=2, threshold=1,
                     buffer_length=1000000):
    '''
    Method to run OSF search in multiple processes simultaneously.
    '''
    if __name__ == '__main__':
        buffer_list = buffer_generator(every_matrix(
            dimension, dimension, full_data), buffer_length)
        pool = Pool(processes=numProcesses)
        identityMat = np.eye(dimension)
        full_data_np = full_data.values
        compiled_chunks = []
        for chunk in buffer_list:
            list_of_combinations = [chunk[i::numProcesses] for i in range(
                numProcesses)]
            result_list = pool.starmap(iterate_RMSs, zip(
                list_of_combinations, repeat(full_data_np.copy()), repeat(
                    identityMat), repeat(threshold)))
            merged_pool = list(chain.from_iterable(result_list))
            compiled_chunks.extend(merged_pool)
        return sorted(compiled_chunks, key=lambda x: x[0])


def main():
    #################################################################
    # setup parser for accepting arguments from the bash shell
    parser = argparse.ArgumentParser(
        description='Multi-process(core) orthogonal set finder.')
    parser.add_argument('-c', '--config_file',
                        help='JSON file containing input configuration.', default=None)
    parser.add_argument('-i', '--input',
                        help='Input file name as .csv. Compounds in columns,'
                        'mutants in rows.')
    parser.add_argument('-o', '--output',
                        help='Output file name as .csv')
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
    parser.add_argument('-b', '--buffer_length',
                        help='Length of list to buffer into memory. Should be'
                        'larger when using more processes. default: 1E6.',
                        default=1000000, type=int)
    parser.add_argument('-e', '--time_testing',
                        help='Use to test speed of algorithm. Restricts printed'
                        ' and csv output. default: off', action='store_true')
    args = vars(parser.parse_args())
    #################################################################
    if args["config_file"] != None:
        conf = json.load(open(args["config_file"]))
    else:
        conf = args
    # run the script printing start and end times and the top five hits at the
    # end.
    if not conf["time_testing"]:
        print('Running a {}x{} matrix search on {} with {} process(es).'
              'Threshold set to {}.'.format(conf["dimension"], conf["dimension"],
                                            conf["input"], conf["processes"],
                                            conf["threshold"]))
    try:
        full_data = pd.read_csv(conf["input"], index_col=0, dtype='float64')
    except FileNotFoundError:
        print("Could not find the file in the specified path: {}. "
              "Please check the file/path.".format(conf["input"]))
        return 0
    except ValueError:
        print("A value in the CSV file was not recognised. Only numbers as floats "
              "or in scientific notation can be processed.")
        return 0
    except IndexError:
        print("It looks like your input matrix is not rectangular. Check for"
              " extraneous commas in the .csv file.")
        return 0
    #    raise
    starttime = datetime.now()
    if not conf["time_testing"]:
        print('Start time: {}'.format(starttime.isoformat()))
    # prep the raw data for processing:
    full_data.index = full_data.index.map(int)  # Allows m numbers to
    # format properly.
    clean_raw_data(full_data)  # Set everything below 1E3 to 1E3.
    # Start the algorithm.
    result = run_multiprocess(full_data, conf["dimension"],
                              numProcesses=conf["processes"],
                              threshold=conf["threshold"],
                              buffer_length=conf["buffer_length"])
    endtime = datetime.now()
    if not conf["time_testing"]:
        print("Done! Top five hits:")
        print(result[:5])
        print("I found {} combinations.".format(str(len(result))))
        print('End time: {}'.format(endtime.isoformat()))
    print('Total calculation time: {}'.format(str(endtime - starttime)))

    if not conf["time_testing"]:
        format_OSF(result, full_data,
                   list_len=conf["length"]).to_csv(conf["output"],
                                                   index=False)  # Write out result.
        print('Result saved to {}'.format(conf["output"]))
    return 1


if __name__ == '__main__':
    exit(main())
