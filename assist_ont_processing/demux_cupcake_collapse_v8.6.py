#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Aggregate isoform counts of similarly grouped isoforms after cupcake collapse to generate final expression matrix
Functions:
    read_input_files(args)
    sum_counts(samples, counts, group, isoform)
    wrapper_sum_counts(samples, counts, group)
Pre-requisites:
    group_txt generated from E.Tseng's collapse_by_isoform.py script (version 8.6)
    expression matrix needs to have matching isoform ID as the collapse group file (under column "union_isoform")
    samples file needs to have matching sample IDs as expression matrix
    

Cupcake collapse takes similar isoforms "collapsed_id" and gives a "new_id"
the representative isoform is the longest isoform 
Aim is to generate an expression matrix with the sum of the counts of the similar isoforms for each collapsed isoform

group: PB.1.1            PB.XX, PB.YY, TALONXX
abundance: isoform        sample_x sample_y
           PB.XX             2       1
           PB.YY             1       1
           TALONXX           1       1
new abundance
           PB.XX             4       3
(PB.XX taken as representive due to most abundant isoform across all the samples)

"""

# packages
import pandas as pd
from gtfparse import read_gtf
import csv
import sys
import argparse
import os
import time
pd.options.mode.chained_assignment = None  # default='warn'


"""
read input files 
:param args: args.sample, args.expression, args.group
:returns sample: dictionary of sample IDs under each dataset
:returns counts: df of raw expression for each sample 
:returns group: df of group file read after cupcake collapse (column1: new_id, column2: collapsed_id)
"""
def read_input_files(args):
    
    # read in file with sample IDs
    samples = pd.read_csv(args.sample)
    
    # read in counts file  
    counts = pd.read_csv(args.expression)
    
    # read in groups file
    group = pd.read_csv(args.cupcake_group, sep = "\t", header = None)
    group.columns = ["new_id","collapsed_id"]
    
    ## check the isoforms in group file are in the abundance file
    # split the group column collapsed_id by commma and combine into one list (iso_lst)
    l = group.apply(lambda row: row["collapsed_id"].split(","), axis=1)
    iso_lst = [item for sublist in l for item in sublist]
    
    # error if isoform in collapsed group has no abundance (i.e. no record in matrix)
    if len(set(iso_lst) - set(counts["union_isoform"])) > 0 :
        print("ERROR: isoforms not in expression matrix!!")
        print(', '.join(set(counts["union_isoform"]) - set(iso_lst) ))
        sys.exit()
        
    #set(abundance["union_isoform"]) - set(iso_lst) 
            
    return samples, counts, group
    

"""
sum the counts for the isoforms that are collapsed across the samples, using the new ID from group file
:param group: read in group file
:param counts: read in expression matrix
:param isoform: isoform  in "new_id" column in group file
:returns results_dict: isoform dictionary of the collapsed counts per sample for the representative most abundant isoform
"""
def sum_counts(samples, counts, group, isoform):
    
    # extract the collapsed isoforms (col_iso) based on the input new id isoform
    col_iso = group[group["new_id"] == isoform]["collapsed_id"].values[0].split(",")
    
    # subset the counts of those collapsd isoforms in the abundance file 
    df = counts[counts["union_isoform"].isin(col_iso)]
    
    # list all the sample IDs (removing any NAs from the input file)
    all_samples = list(samples[samples.columns[0]].values) + list(samples[samples.columns[1]].values)
    all_samples = [s for s in all_samples if str(s) != 'nan']
    
    # sum the counts for each column (i.e. each sample) after subsetting
    dfsum = pd.DataFrame(df[["union_isoform"] + all_samples].sum(axis=0))
    dfsum.columns = ["count"]

    # create a dictionary
    results_dict = dfsum['count'].to_dict()
    
    ## determine the most abundant isoform
    # sum across all the samples for each isoform
    df["sumFL"] = df[all_samples].sum(axis=1)
    abundant_iso = df.loc[df['sumFL'].idxmax()]["union_isoform"]
    
    # replace union isoform with the most abundant isoform 
    # note this would be the concatenated collapsed isoforms at this point
    results_dict.update({"union_isoform" : abundant_iso})

    return results_dict


"""
apply the sum_counts for each isoform to the all the isoforms in the group file
:param group:  read in group file
:returns final: final expression matrix of collapsed counts across isoforms
"""
def wrapper_sum_counts(samples, counts, group):
    # apply sum_counts() to each row in the group file
    output = group.apply(lambda row: sum_counts(samples, counts, group, row["new_id"]), axis=1)
    
    # append each dictionary from each isoform (output of sum_counts) to a list
    # convert list to dataframe
    final = []
    for i in output: 
        final.append(i)
    final = pd.DataFrame(final)
    final.set_index('union_isoform')
  
    return(final)


def main():
    parser = argparse.ArgumentParser(description="Identifying Transcripts from ONT and Iso-Seq Targeted Transcriptome for annotation")
    parser.add_argument("expression", help='\t\traw expression .txt file.')
    parser.add_argument("cupcake_group", help='\t\tgroup.txt file generated from cupcake collapse')
    parser.add_argument('sample', help='\t\tTxt file containing sample IDs; column 1 - Iso-Seq, column 2 - ONT')
    parser.add_argument('--out', help='\t\tPath and name of output file')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output abundance file; Default: <group_prefix>_fl_count.csv.')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output gtf. Default: group.txt.')
   
    args = parser.parse_args()
    
    if args.dir is None:
      args.dir =  os.path.dirname(args.cupcake_group)
      
    if args.output is None:
      args.output = os.path.basename(args.cupcake_group.replace(".group.txt",""))
    
    print("Reading input files...")
    samples, counts, group = read_input_files(args)
    
    print("Aggregating reads of collapsed isoforms across samples...")
    final = wrapper_sum_counts(samples, counts, group)
    
    filename = args.dir + "/" + args.output + "_fl_count.csv"
    print("Writing output to", filename)
    final.to_csv(filename, index = False) 
    
    print("All Done")
    
if __name__ == "__main__":
    main()
