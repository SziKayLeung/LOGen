#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: identify raw ont reads collapsed using TALON
Input: 
    talon_annot_transcript_id of interest for subsetting 
    talon_annot_read.tsv generated from TALON
    ont_gtf of combined aligned raw reads 
Output:
    subsetted gtf of aligned raw reads which were clustered to talon_annot_transcript_id of interest
Functions:
    subset_gtf_with_raw_reads
    identify_raw_reads
Pre-requisite:
    original raw read id in talon_annot_read.tsv is documented in column 1 
    talon_annot_transcript_id in talon_annot_read.tsv is documented in column 13 

---- Input:    
==> talon_annot_transcript_id
XXXX ....... TALONT0001234
AAAA ....... TALONT0001234
YYYY ....... TALONT0001111

==> ont_gtf
chr ...... transcript_id "XXXX"
chr ...... transcript_id "AAAA"
chr ...... transcript_id "YYYY"

==> isoform of interest = TALONT0001234

---- Output:
==> TALONT0001234_rawread.gtf 
chr ...... transcript_id "XXXX"
chr ...... transcript_id "AAAA"


"""

# packages
import pandas as pd
from gtfparse import read_gtf
import csv
import argparse
import os


"""
subset_gtf_with_raw_reads
write line into write_df if the strings in raw_reads matches the input_gtf
:params args
:params raw_reads = list of strings for matching i.e. [XXXX AAAA]
:params write_df = opened output file to write lines to
""" 
def subset_gtf_with_raw_reads(args, raw_reads, write_df):
    # open gtf and loop through each line
    with open(args.gtf) as file:
        for line in file:
            if any(read in line for read in raw_reads):
                write_df.write(line)
                

"""
identify_raw_reads
search the isoform of interest in the talon_annot.tsv and record the original read id
grep the original read id from gtf and output
:params args
:writes <isoform>_rawread.gtf
""" 
def identify_raw_reads(args):
    print("Extracting ONT raw reads collapsed using TALON...")
    
    # loop through each isoform of interest to output
    for i in args.isoform:
        print("Collapsed to", i)
        
        raw_reads = []
        with open(args.talon_annot) as annot_file:
            for line in annot_file:
                # identify the original raw read and transcript_id in annot_reads.tsv 
                # 1st and 13th column
                raw_read = line.split("\t")[0]
                transcript_id = line.split("\t")[14]
                if i == transcript_id:
                  print(transcript_id)
                  raw_reads.append(raw_read)
        
        # subset gtf with raw_reads list
        filename = args.dir + "/" + i +"_rawread.gtf"
        print("Writing to", filename)
        gtf_out_df = open(filename,'w')
        subset_gtf_with_raw_reads(args, raw_reads, gtf_out_df)
        gtf_out_df.close()
    
  
def main():
    parser = argparse.ArgumentParser(description="Identify raw reads collapsed using TALON")
    parser.add_argument('talon_annot', help='\t\tTALON output annot_reads.tsv.')
    parser.add_argument('gtf', help='\t\tInput gtf of (combined) ONT reads')
    parser.add_argument('--isoform', nargs='+',required=True, help='\t\tList of TALON collapsed isoform ID to extract; isoform <talon_id1> <talon_id2>')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output gtf. Default: Directory of annot_reads.tsv.')

    args = parser.parse_args()
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.talon_annot)
        
    identify_raw_reads(args)
    print("All Done")

if __name__ == "__main__":
  main()
