#!/usr/bin/env python3

import os
import pandas as pd
import argparse
from functools import reduce

def read_demux(args):
    # List all files with a specific pattern
    demux_names_files = [f for f in os.listdir(args.input_directory) if f.endswith(args.chr + "_fl_count.csv")]
    print("Total number of files to read:", len(demux_names_files))
    
    # Load and process each CSV file
    demux_files = []
    for file_name in demux_names_files:
        file_path = os.path.join(args.input_directory, file_name)
        print("***** Reading:", file_path)
        df = pd.read_csv(file_path)
        
        # Convert `id` column to string and set it as row index
        df['id'] = df['id'].astype(str)
        df.set_index('id', inplace=True)
        
        demux_files.append(df)

    return(demux_files)

def merge_demux(demux_files):
    
    # Merge the DataFrames using a full outer join
    print("***** Merging dataframe")
    merged_df = reduce(lambda left, right: pd.merge(left, right, on='id', how='outer'), demux_files)
    merged_df.fillna(0, inplace=True)
    merged_df = merged_df.astype(int)

    return(merged_df)

def sum_nreads(merged_df):
    # sum number of reads and samples across all datasets
    nreads = merged_df.sum(axis=1)
    nsamples = (merged_df != 0).sum(axis=1)
    # sum number of reads and samples across targeted datasets
    targeted_df = merged_df.filter(like='Targeted')
    targeted_nreads = targeted_df.sum(axis=1)
    targeted_nsamples = (targeted_df != 0).sum(axis=1)
    # sum number of reads and samples across whole datasets   
    whole_df = merged_df.filter(like='Whole')
    whole_nreads = whole_df.sum(axis=1)
    whole_nsamples = (whole_df != 0).sum(axis=1)
    # populate info
    merged_df["nreads"] = nreads
    merged_df["nsamples"] = nsamples
    merged_df["targeted_nreads"] = targeted_nreads
    merged_df["targeted_nsamples"] = targeted_nsamples
    merged_df["whole_nreads"] = whole_nreads
    merged_df["whole_nsamples"] = whole_nsamples

    return(merged_df)

def main():
    parser = argparse.ArgumentParser(description="Extracting polyA and polyT sequences")
    parser.add_argument('-i', "--input_directory", help='\t\tInput directory containing demux files.')
    parser.add_argument("--chr", help='\t\tChromosome to read in')
    parser.add_argument('-o_name', "--output_name", help='\t\tOutput name')
    parser.add_argument('-o_dir', "--output_dir", help='\t\tOutput directory')

    args = parser.parse_args()

    demux_files = read_demux(args)
    merged_df = merge_demux(demux_files)
    merged_df = sum_nreads(merged_df)

    # output
    output_name = args.output_dir + "/" + args.output_name + "_" + args.chr + "all_fl_count.csv"
    print("***** Writing output to:", output_name)
    merged_df.to_csv(output_name)


if __name__ == "__main__":
    main()
    
    

