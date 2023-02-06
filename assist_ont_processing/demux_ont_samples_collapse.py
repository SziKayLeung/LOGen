#!/usr/bin/env python3

import pandas as pd
import csv
import argparse
import sys


def main(args):

    # read in files
    print("Reading in files...")
    readstat = pd.read_csv(args.stat, sep = "\t")
    samples = pd.read_csv(args.read)
    
    # merge read_stat and samples files using "id" column
    print("Merging by id...")
    merged = pd.merge(readstat, samples, how = "outer")
    
    # replace pb.id NA values with 0 
    # example where there is an entry in the readstat file but not in the samples file
    mergedna = merged[merged['pbid'].isna()]
    mergedna.to_csv(args.o + "/missingID_demux_fl_count.txt", index = False)
    
    merged['pbid'] = merged['pbid'].fillna(0)
    
    
    # tally and pivot 
    merged_tally = merged.groupby(['pbid','primer']).count()
    merged_tally = merged_tally.reset_index()
    merged_tally_pivot = merged_tally.pivot(index="pbid", columns="primer", values="id")
    
    if args.sample is not None:
      print("Replace batch and barcode names with sample names")
      sample = pd.read_csv(args.sample, sep = "\t")
      merged_tally_pivot.rename(columns=sample.set_index('BatchBC')['#Sample'], inplace=True)
    
    # replace NA with 0
    # reset index and rename column
    merged_tally_pivot.fillna(0, inplace=True)
    merged_tally_pivot.reset_index(inplace=True)
    merged_tally_pivot.rename({'pbid': 'isoform'}, axis=1, inplace=True)
  
    # write output
    merged_tally_pivot.to_csv(args.o + "/demux_fl_count.csv", index = False)
    
    print("All Done")
    

if __name__ == "__main__":
  
    parser = argparse.ArgumentParser(description="Creating a cluster_report.csv format as input for cupcake tofu get_abundance_post_collapse.py")
    parser.add_argument("--stat", help='\t\tread_stat.txt generated from Cupcake/Isoseq3 collapse')
    parser.add_argument('--read', help='\t\tID file of reads')
    parser.add_argument('--o', help='\t\toutput_dir')
    parser.add_argument('--sample', required=False, default=None, help='\t\ttsv file of samples to replace <#sample> <BatchBC> (optional)')
  
    args = parser.parse_args()
  
    main(args)
