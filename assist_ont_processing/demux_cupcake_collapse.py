#!/usr/bin/env python3


"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Date: Jan 2023
Aim: Uses read_stat.txt from IsoSeq3-collapse and unique ONT read ID to generate an abundance file at sample-level

---- Input:
==> <sample>.read_stat.txt = csv file documenting original reads clustered
id,pbid
AA,PB.1.1
BB,PB.1.1
transcript/1,PB.1.1

arg.dataset = ont
==> <all_sample>.sample_id.csv = csv file documenting origin of reads for all samples
id,primer
AA,Sample1
BB,Sample2

arg.dataset = isoseq
==> <all_sample>.sample_id.csv = csv file documenting origin of reads for all samples
cluster_id,id,primer
transcript/1,mXX/ccs,Sample3
transcript/1,mYY/ccs,Sample4

---- Output:
arg.dataset = ont
==> demux_fl_count.csv
        Sample1  Sample2
PB.1.1    1        1

arg.dataset = isoseq
        Sample3  Sample4
PB.1.1    1        1

"""

# packages
import pandas as pd
import os
import csv
import argparse
import sys

"""
demux_ont
:params args; args.read_stat & args.sample_id already read in
:writes demux_ignored_id.txt
:writes demux_fl_count.csv
:return merged_tally_pivot = df of demux_fl_count.csv
""" 
def demux(args):
    '''
    1. Merge files using "id" column or "cluster_id" and "id" i.e. read ID; keeping only common reads documented in both files
    for ONT dataset, merge using common column "id"
    NB all reads in sample_id.csv will be read_stat.txt, as sample_id.csv is derived from the .fa that should be used for IsoSeq3 collapse
    
    for Iso-Seq dataset, merge using id and cluster_id
    Iso-Seq datasets, unlike ONT datasets, processed with Iso-Seq3 cluster 
    Iso-Seq3 cluster labels ID as "transcript/X"
    Iso-Seq collapse with Iso-Seq datasets therefore uses Iso-Seq3 cluster ID (rather than the original ccs ID)
    NB keep only reads in read_stat file from collapse
    
    2. Writes reads that have been discarded by Iso-Seq3 collapse 
    3. Tally the number of reads across isoform and sample
    4. Replace the sample column if --sample argument provided
    '''
    
    # merge read_stat and samples files; columns to merge dependent on dataset
    if args.dataset == "ont":
        print("Processing ONT dataset")
        merged = pd.merge(args.read_stat, args.sample_id, how = "outer", on = "id")
    else:
        print("Processing Iso-Seq dataset")
        merged = pd.merge(args.read_stat, args.sample_id[["cluster_id","primer"]], how = "outer", left_on = "id", right_on = "cluster_id")
    
    # document reads that are discarded by Iso-Seq3 collapse
    # i.e. entry in sample_id.csv but not in read_stat.txt
    mergedna = merged[merged['pbid'].isna()]
    
    # replace pb.id NA values with 0 
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
    # reset index and rename column for input to SQANTI
    # drop 0 index as this represents reads that were discarded during collapse
    merged_tally_pivot.fillna(0, inplace=True)
    merged_tally_pivot.reset_index(inplace=True)
    merged_tally_pivot.rename({'pbid': 'id'}, axis=1, inplace=True)
    merged_tally_pivot.drop([0], inplace = True)
  
    # write output
    mergedna.to_csv(args.dir + "/" + args.output + "_ignored_id.txt", index = False)
    output_count = args.dir + "/" + args.output + "_fl_count.csv"
    print("Writing output files to:", output_count)
    merged_tally_pivot.to_csv(output_count, index = False)
    
    return(merged_tally_pivot)


def main():
    parser = argparse.ArgumentParser(description="ONT demultiplex and retrieve read counts from IsoSeq3 collapse output")
    parser.add_argument("read_stat", help='\t\tread_stat.txt generated from Cupcake/Isoseq3 collapse')
    parser.add_argument('sample_id', help='\t\tcsv file of reads and sample; <id,primer> for ont\t <cluster_id,id,primer> for isoseq')
    parser.add_argument('--dataset', default = "ont", choices=['ont', 'isoseq'], help='\t\tdataset:"ont" or "isoseq"')
    parser.add_argument('--sampleLevel', default = True, help='\t\tmerging at a sample level i.e headers = <read_id>"')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output abundance file; Default: demux_fl_count.csv.')
    parser.add_argument('-s','--sample', required=False, default=None, help='\t\ttsv file of samples to replace <#sample> <BatchBC> (optional)')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of read_stat.txt.')
  
    args = parser.parse_args()
    if args.output is None:
        args.output = "demux" 
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.read_stat)
        
    # read in files
    print("Reading in read stat file...")
    args.read_stat = pd.read_csv(args.read_stat)
    
    print("Reading in sample id file...")
    _, extension = os.path.splitext(args.sample_id)
    if ".txt" == extension:
        print("as txt file")
        args.sample_id = pd.read_csv(args.sample_id, sep = "\t")
    else:
        print("as csv file")
        args.sample_id = pd.read_csv(args.sample_id)

    demux(args)
    print("All Done")
    

if __name__ == "__main__":
    main()
