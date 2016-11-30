'''
Runs OSF search in the terminal. Requires the orthogonal_set_finder.py file.
'''
from orthogonal_set_finder import *
import argparse
from datetime import datetime
import pandas as pd
from multiprocessing import Pool
from itertools import repeat, chain

def run_multiprocess(full_data, dimension, numProcesses=2):
    '''
    Method to run OSF search in multiple processes simultaneously.
    '''
    if __name__ == '__main__':
        # [ lst[i::n] for i in xrange(n) ] Splits up lst into n segments
        list_of_combinations = [ list(every_matrix(dimension,dimension,full_data))[i::numProcesses] for i in range(numProcesses) ]
        pool = Pool(processes=numProcesses)
        result_list = pool.starmap(iterate_RMSs, zip(list_of_combinations, repeat(full_data.copy())))
        merged = list(chain.from_iterable(result_list))
        return sorted(merged, key=lambda x: x[0])

#################################################################
# setup parser for accepting arguments from the bash shell
parser = argparse.ArgumentParser(description='Multi-process(core) orthogonal set finder.')
parser.add_argument('-i','--input', help='Input file name as .csv. Compounds in columns, mutants in rows.', required=True)
parser.add_argument('-o','--output', help='Output file name as .csv', required=True)
parser.add_argument('-d','--dimension', help='Dimension of search. Default: 2.', default=2, type=int)
parser.add_argument('-p','--processes', help='Number of processes to spawn. Default: 1.', default=1, type=int)
parser.add_argument('-l','--length', help='Length of result list to print to csv. default: 1000.', default=1000, type=int)
args = parser.parse_args()

# run the script printing start and end times and the top five hits at the end.
print('Running a {}x{} matrix search on {} with {} process(es).'.format(args.dimension, args.dimension, args.input, args.processes))
try:
    full_data = pd.read_csv(args.input, index_col=0, dtype='float64')
except:
    print("Something went wrong with the import of {}. Please check the file/path.".format(args.input))
    raise
starttime = datetime.now()
print('Start time: {}'.format(starttime.isoformat()))
clean_raw_data(full_data) # Set everything below 1E3 to 1E3.
result = run_multiprocess(full_data, args.dimension, numProcesses = args.processes)
print("Done! Top five hits:")
print(result[:5])
endtime = datetime.now()
print('End time: {}'.format(endtime.isoformat()))
print('Total calculation time: {}'.format(str(endtime - starttime)))
# May have to tweak the .to_csv() arguments to get rid of float formatting of
# mutant numbers and extra quotes that appear around compounds.
format_OSF(result, full_data, list_len=args.length).to_csv(args.output, index=False) # Write out result.
print('Result saved to {}'.format(args.output))
#################################################################
