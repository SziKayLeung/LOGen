#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Date: Jan 2023
Aim: Uses read_stat.txt from merged Iso-Seq and ONT dataset to generate an abundance file at sample-level across both datasets

Premise:
Merge ONT dataset (fasta generated from Transcript Clean) + Iso-Sea dataset (fasta generated from Iso-Seq3 cluster)
Able to differentiate ONT and Iso-Seq reads, as Iso-Seq reads have "Transcript" in read prefix 
1. Therefore split the read_stat.txt generated from Iso-Seq Collapse into isoforms that have Iso-Seq reads collapsed, 
and isoforms with only ONT reads
2. run the demux() from demux_cupcake_collapse.py on ONT and Iso-Seq read_stat, supplying respective sample_id.csv files
3. merge the two abundance files generated by the isoform ID

"""

# packages
from Bio import SeqIO
import os
import glob
import argparse
from argparse import Namespace
import pandas as pd
import sys

# custom script
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sys.path.append(LOGEN + "compare_datasets")
import demux_cupcake_collapse as dcc
import dataset_identifer as di


"""
split_reads_by_datset
split read_stat.txt file into two read_stat.txt consisting of isoforms with Iso-Seq and ONT reads respectively
:params args
:return isoseq_read_stat = read_stat file of Iso-Seq-derived reads
:return ont_read_stat = read_stat file of ONT-derived reads
""" 
def split_reads_by_datset(args):
    '''
    merged_read_stat file generated from Iso-Seq collapse of ONT and Iso-Seq reads (clustered)
    id            pbid
    XXX           PB.1.1
    YYY           PB.1.1
    transcript/1  PB.1.1
    
    Iso-Seq clustered reads all have prefix "transcript/<number>"
    Differentiate ONT and Iso-Seq reads by "transcript" in readID
    
    ==> isoseq_read_stat
    transcript/1 PB.1.1
    ==> ont_read_stat
    XXX          PB.1.1
    YYY          PB.1.1
    '''
    merged_read_stat = pd.read_csv(args.merged_read_stat, sep = "\t")
    
    # subset Iso-Seq reads from read_stat file
    isoseq_read_stat = merged_read_stat[merged_read_stat.id.str.contains('transcript')]
    
    # subset ONT reads as not containing "transcript"
    ont_read_stat = merged_read_stat[~merged_read_stat.id.str.contains('transcript')]
    
    return isoseq_read_stat, ont_read_stat


"""
demux_dataset
generate abundance file of read counts from Iso-Seq and ONT dataset
:params args
:writes demux_ont_fl_count.csv = isoform abundance from ONT reads by sample
:writes demux_iso_fl_count.csv = isoform abundance from Iso-seq reads by sample
:writes demux_fl_count.csv = merged isoform abundance of Iso-Seq and ONT reads by sample
:return merged_abundance = df of demux_fl_count.csv
""" 
def demux_dataset(args):
    '''
    ---- Input:
    ==> ont_sample_id
    id,primer
    XXX, S1
    YYY, S2
    
    ---- Output:
    ==> ont_abundance
              S1    S2
    PB.1.1    1      1
    ==> merged_abundance
             ont_S1  ont_S2
    PB.1.1 
    '''
    
    # split reads by dataset
    isoseq_read_stat, ont_read_stat = split_reads_by_datset(args)
    
    # read sample ID files
    ont_sample_id = pd.read_csv(args.ont_sample_id)
    iso_sample_id = pd.read_csv(args.iso_sample_id)
    
    # arguments to demultiplex ONT dataset
    ont_demux_args = Namespace(read_stat = ont_read_stat, sample_id = ont_sample_id, dataset = "ont", 
                               sample = args.ont_sample, dir = args.dir, output = args.output + "_ont")
    
    iso_demux_args = Namespace(read_stat = isoseq_read_stat, sample_id = iso_sample_id, dataset = "isoseq", 
                               sample = args.iso_sample, dir = args.dir, output = args.output + "_iso")

    # run demux() in demux_ont_samples_collapse.py 
    print("Demultiplexing....")
    ont_abundance = dcc.demux(ont_demux_args)
    iso_abundance = dcc.demux(iso_demux_args)
    
    # add dataset names to columns in abundance files for later merging
    ont_abundance.columns = ["isoform"] + ["ONT_" + i for i in ont_abundance.columns[1:]]
    iso_abundance.columns = ["isoform"] + ["Iso-Seq_" + i for i in iso_abundance.columns[1:]]
    
    # merge ont and iso-seq abundance
    # keep all isoforms in both datasets
    print("Merging abundance....")
    merged_abundance = pd.merge(ont_abundance,iso_abundance,on="isoform",how="outer")
    
    # create column of sum of read count in each dataset
    merged_abundance["ONT_sum_FL"] = merged_abundance[merged_abundance.filter(like='ONT').columns].sum(axis=1)
    merged_abundance["Iso-Seq_sum_FL"] = merged_abundance[merged_abundance.filter(like='Iso-Seq').columns].sum(axis=1)
    
    # replace NA in columns (where unique isoforms in dataset) with 0
    merged_abundance.fillna(0, inplace=True)
    
    # remove isoforms labelled as "0"
    # demux() in demux_ont_samples_collapse.py tallies all isoforms that are discarded in IsoSeq collapse
    # count is given under isoform ID "0"
    merged_abundance = merged_abundance.loc[merged_abundance["isoform"]!=0,]
    
    # set "isoform" to index to create array
    merged_abundance.set_index("isoform",inplace=True)
    
    # convert array of floats to integers, given count data
    merged_abundance = merged_abundance.astype(int)
    
    # classify dataset based on count data    
    merged_abundance["dataset"] = merged_abundance.apply(
        lambda row: di.classifier_dataset_bycount(row['ONT_sum_FL'], row['Iso-Seq_sum_FL'],"ONT","Iso-Seq"),
        axis=1)
    
    # write output
    final_output = args.dir + "/" + args.output + "_fl_count.csv"
    print("Writing output files to:", final_output)
    merged_abundance.to_csv(final_output)
    
    return(merged_abundance)


def main():
    parser = argparse.ArgumentParser(description="Demultiplex merged cluster_report.csv using sample flnc read id")
    parser.add_argument("merged_read_stat", help='\t\tread_stat.txt from collapsing isoseq + ont datasets')
    parser.add_argument('ont_sample_id', help='\t\ttxt file of ONT reads across samples; <id,primer>')
    parser.add_argument('iso_sample_id', help='\t\ttxt file of Iso-Seq reads across samples; <cluster_id,id,primer>')
    parser.add_argument('--ont_sample', required=False, default=None, help='\t\ttsv file of ONT samples to replace <#sample> <BatchBC> (optional)')
    parser.add_argument('--iso_sample', required=False, default=None, help='\t\ttsv file of Iso-Seq samples to replace <#sample> <BatchBC> (optional)')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for demuxed_cluster_report.csv; Default: demux')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of merged read_stat.txt.')
  
    args = parser.parse_args()
    
    if args.output is None:
        args.output = "demux" 
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.merged_read_stat)
        
    demux_dataset(args)
    print("All Done")
    

if __name__ == "__main__":
    main()
