#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: identify raw reads collapsed using Iso-Seq3 collapse/cupcake collapse
Input: 
    --isoform of interest for subsetting 
    group.txt generated from Iso-Seq3 collapse/cupcake collapse
    gtf of reads; equivalent of .fasta input for collapse
Output:
    subsetted gtf of aligned raw reads which were collapsed to isoform of interest
    gtf of representiative isoform of interest (--rep)
Functions:
    read_input
    identify_raw_reads
Pre-requisite:
    Based on group.txt generated from cupcake collapse (see below for format)

---- Input:
==> group.txt
PB.1.1 rawread1,rawread2,rawread3
PB.1.2 rawread4

"""

# packages
import pandas as pd
from gtfparse import read_gtf
import argparse
import sys
import os

import warnings
warnings.filterwarnings("ignore")

# LOGen modules
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sys.path.append(LOGEN + "miscellaneous")
import read_write as rw


"""
read_input
read in group.txt file from cupcake collapse, and original gtf used for cupcake collapse
grep the original read id from gtf and output
:params args
:return cols_group = read in cupcake collapse group file
:return input_gtf = read in original aligned gtf
""" 
def read_input(args):
  
    # read in gtf file used for collapsing cupcake
    print("**** Reading in:", args.gtf)
    input_gtf = read_gtf(args.gtf)
    
    # read in cupcake output files
    print("**** Reading in:", args.cupcake_group)
    col_group = pd.read_csv(args.cupcake_group, delimiter='\t', header = None)
    col_group.columns = ["final_isoform", "iso_collapsed"]
    
    return col_group, input_gtf
  

"""
identify_raw_reads
subset the input gtf by isoform of interest after extracting the raw reads from group.txt
:params args
:writes <isoform>_pre_collapsed.gtf = gtf of raw reads collapsed into isoform of interest
:writes <isoform>_rep_collapsed.gtf = gtf of representative isoform from collapsing
:return subset_input_gtf = read in original aligned gtf
"""
def identify_raw_reads(args):
    # read input
    col_group, input_gtf = read_input(args)
    
    print("Extracting raw reads collapsed using Iso-Seq collapse/cupcake collapse...")
    
    for i in args.isoform:
        print("Collapsed to", i)
        
        # subset the 2nd column of group.txt file, string contains to work with cupcake collapse output
        df = col_group[col_group["iso_collapsed"].str.contains(i)]
        
        # extract all the collapsed original isoforms into one list
        lst = []
        for iso in df["iso_collapsed"]: lst.append(str(iso))
        
        # split the list by comma; given <isoform1,isoform2,isoform3>
        # and subset input gtf
        subset_input_gtf = input_gtf[input_gtf['transcript_id'].isin(str(lst[0]).split(','))]
        output_gtf_file = args.dir + "/" + i + "_pre_collapsed.gtf"
        print("Writing to", output_gtf_file)
        rw.writeGTF(subset_input_gtf, output_gtf_file)
        
        # also extract the representative gtf if cupcake_gtf provided
        if args.cupcake_gtf is not None:
            subset_collapsed_gtf = args.cupcake_gtf[args.cupcake_gtf["transcript_id"] == df["final_isoform"].values[0]]
            output_rep_file = args.dir + "/" + i + "_rep_collapsed.gtf"
            rw.writeGTF(subset_collapsed_gtf, output_rep_file)
        
    return subset_input_gtf
  
  
def main():
    parser = argparse.ArgumentParser(description="Identify raw reads collapsed using Cupcake collapse")
    parser.add_argument('cupcake_group', help='\t\tIso-Seq collapse/cupcake collapse group.txt.')
    parser.add_argument('gtf', help='\t\t Input gtf corresponding to input fasta for cupcake collapse.')
    parser.add_argument('--rep', default=False, action="store_true", help='\t\tGenerate represented isoform from cupcake collapse')
    parser.add_argument('--isoform', nargs='+',required=True, help='\t\tList of isoform ID to extract; isoform <isoform_id1> <isoform_id2>')
    parser.add_argument('--cupcake_gtf', required=False, help='\t\tIso-Seq collapse/cupcake collapse group.gff (optional).')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output gtf. Default: group.txt.')
    
    args = parser.parse_args()
    print("************ Extract original sequences...", file=sys.stdout)
    
    if args.dir is None:
      args.dir =  os.path.dirname(args.cupcake_group)
    
    # generate representative transcript by using collapsed.gtf
    # assume cupcake gtf is in collapsed folder
    if args.rep:
      collapsed_gtf = args.cupcake_group.replace(".group.txt",".gff")
      print("Extract representative isoform from collapsing", collapsed_gtf)
      args.cupcake_gtf = read_gtf(collapsed_gtf)
        
    identify_raw_reads(args)
    print("**** All Done! ****")
  
  
if __name__ == "__main__":
    main()
