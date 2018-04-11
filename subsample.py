'''
Accepts csv number, and fractions as arguments and writes random subsamples out.
'''

import argparse
import json
import random
import pandas as pd

def random_df_sample(df, fraction):
    indicies = list(df.index)
    sample = random.sample(indicies, int(fraction * len(indicies)))
    return df.loc[sample]

def main():
    #################################################################
    # setup parser for accepting arguments from the bash shell
    parser = argparse.ArgumentParser(
        description='Subsampler')
    parser.add_argument('-i', '--input',
                        help='Input file name as .csv. Compounds in columns,'
                        'mutants in rows.')
    parser.add_argument('-o', '--output',
                        help='Output file name as .csv')
    parser.add_argument('-f', '--fraction',
                        help='Fraction of subsample.', default=.2, type=float)
    parser.add_argument('-n', '--file_number',
                        help='Number of files to output.',
                        default=5, type=int)
    args = vars(parser.parse_args())
    #################################################################
    df = pd.read_csv(args["input"], index_col=0)
    sub_dfs = [random_df_sample(df, args["fraction"]) for _ in range(args["file_number"])]
    for i, sub_df in enumerate(sub_dfs):
        sub_df.to_csv(args["output"][:-4] + str(i+1) + '.csv')
    return 1


if __name__ == '__main__':
    exit(main())
