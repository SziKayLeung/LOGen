#!/usr/bin/env python

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Subset gtf based on transcript ID and read output
Functions:
    subset_gtf
"""

# packages
import os 
import sys
import glob
import argparse
import pandas as pd
from Bio import SeqIO
from itertools import chain


"""
Subset fasta based sequence match
:param seq_search: string of sequence to match  
:param input_fasta: path of fasta for subsetting
:param name: string of name to be included in the output
:writes subsetted fasta if matches
:returns allcounts, counts = number of total reads, number of reads captured
"""      
def subset_fasta_byseq(seq_search, input_fasta, output_dir, name):
    
    # output names 
    basename = '.'.join(os.path.basename(input_fasta).split(".")[:-1])
    outputfile = output_dir + "/" + basename + "_" + name +".fa"
    
    # read ID
    # number of total reads
    # number of reads captured
    readID = []
    allcounts = int()
    counts = int()
    
    # loop through sequence and find match
    fasta_sequences = SeqIO.parse(open(input_fasta),'fasta')
    with open(outputfile, "w") as f:
        for record in fasta_sequences:
            allcounts = allcounts + 1
            if seq_search in record.seq:
                SeqIO.write([record], f, "fasta")
                readID.append(record.description)
                counts = counts + 1
    
    # delete empty fasta file if not matches
    if os.path.getsize(outputfile) == 0:
        os.remove(outputfile)
    
    print("Writing to:", basename + "_" + name +".fa")
                
    return readID, allcounts, counts


def parse_via_fasta(args):
    readOutput = {}
    countOutput = {}
    for file in glob.glob(args.fasta + "/" + "*" + args.index + "*"):
        print("Reading in", file)
        samplename = os.path.basename(file).split(".")[0]
        readID, allcounts, counts = subset_fasta_byseq(args.seq,file,args.dir,args.output)
        for r in readID:
          readOutput[r] = samplename
        countOutput[samplename] = [counts, allcounts]
    
    finalcounts = pd.DataFrame.from_dict(countOutput, orient='index', columns = ['matched_reads','all_reads'])
    finalcounts.to_csv(args.dir + "/" + args.output + "_all_counts.csv", index_label="sample")
    
    finalreads = pd.DataFrame.from_dict(readOutput, orient='index', columns = ['sample'])
    finalreads.to_csv(args.dir + "/" + args.output + "_all_reads.csv", index_label="matched_reads")

    
def main():
    parser = argparse.ArgumentParser(description="Subset fasta and gtf from list of ID")
    parser.add_argument('--fasta', help='\t\tInput directory with fasta files')
    parser.add_argument('--seq', help='\t\tString of sequence to search')
    parser.add_argument('-i','--index', required=False, default = ".fa", help='\t\tSpecific string in fasta files to capture in input directory')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output file; Default: out.')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: directory of input file.')
  
    args = parser.parse_args()
    
    if args.output is None:
        args.output = "out" 
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.fasta)
    print("Writing output files to:", args.dir)
    
    parse_via_fasta(args)
    
    print("All Done")
    

if __name__ == "__main__":
    main()
