#!/usr/bin/env python3


"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Date: Jan 2023
Aim: generate supplement files for ont pipeline
1. <sample>.cluster_report.csv (equivalent format to Iso-Seq pipeline) to run get_abundance_post_collapse.py
2. <prefix>_sample_id.csv for demultplexing

Example of input and output files generated:
---- Input:
.../Batch1/S1.fasta
> A
TGAT....
> B
AGCA....

.../Batch1/S2.fasta
> C
TGAT....

---- Output:
.../Batch1/S1.cluster_report.csv 
cluster_id,read_id,read_type
transcript/A,A,FL
transcript/B,B,FL

.../Batch1/S2.cluster_report.csv 
cluster_id,read_id,read_type
transcript/C,C,FL

.../Batch1/<all_samples>_sample_id.csv
A,S1
B,S1
C,S2

"""

# packages
import argparse
import os
import glob
from Bio import SeqIO


"""
list_named_fasta
to return specified sample.fa, sample name must in the 1st part of prefix (i.e S1_xxxx.fasta)
:returns list: fasta files in input directory 
"""
def list_named_fasta(args):
    # Store the names of fasta files
    fasta_files = []
    
    # Use recursive glob to search for files in all subdirectories
    search_pattern = args.fasta + "/**/*" + args.index + "*"
    
    # Find all files that match the pattern recursively
    for file in glob.glob(search_pattern, recursive=True):
        
        # If the selected samples are specified
        if args.samples is not None:
            for sample in args.samples:
                if sample in file:
                    print("Including:", file, "\n")
                    fasta_files.append(file)
        else:
            #print("Including:", file, "\n")
            fasta_files.append(file)
    
    return fasta_files


"""
extract_reads_output_cluster
:params args
:params infasta = str of input fasta file
:params outsampledf = open <all_samples>_sample_id.csv to append read id and sample
:writes <sample>_cluster_report.csv
:writes <all_samples>_sample_id.csv 
""" 
def extract_reads_output_cluster(args, infasta, outsampledf):
    '''
    infasta = ../../Batch1/S1_input.fasta
    
    1. records read id of each infasta, and writes to:
    <path_to_infasta>/<prefix_of_infasta>_cluster_report.csv if args.batch is off or
    <path_to_infasta>/<dirname_of_infasta>_<prefix_of_infasta>_cluster_report.csv if args.batch is on
    i.e. 
    args.batch = off: out_cluster = ../../Batch1/S1_cluster_report.csv
    args.batch = on: out_cluster = ../../Batch1/Batch1S1_cluster_report.csv

    2. appends read id of each infasta to outsampledf (open file created outside of function)
    i.e. 
    inputfasta = ../../Batch1/S1_input.fasta
    args.batch = off:
        A,S1
        B,S1
    args.batch = on
        A,Batch1S1
        B,Batch1S1
    '''
    
    print("Processing:", infasta)
    fasta_sequences = SeqIO.parse(open(infasta),'fasta')
    
    # determine whether to include previous directory name into output
    if args.batch:
        print("Samples are processed in batches, including directory name in output")
        infastadir = os.path.dirname(infasta)
        batchname = os.path.basename(infastadir)
        sample = batchname + os.path.basename(infasta).split(".")[0]
    else:
        sample = os.path.basename(infasta).split(".")[0]
    
    # create <sample>_cluster_report.csv
    out_cluster_name = args.dir + "/" + sample + "_cluster_report.csv"
    print("Writing to:", out_cluster_name,"\n")
    out_cluster = open(out_cluster_name,"w")
    out_cluster.write("cluster_id,read_id,read_type\n")
    
    out_sample_name = args.dir + "/" + sample + "_sample_id.csv"
    out_sample = open(out_sample_name, "w")
    out_sample.write("id,primer\n")
    
    # write each fasta ID in fasta sequences
    for fasta in fasta_sequences: 
        out_cluster.write("transcript/" + fasta.id + "," + fasta.id + ",FL\n")
        out_sample.write(str(fasta.id) + "," + sample + "\n")
        outsampledf.write(str(fasta.id) + "," + sample + "\n")
    
    # close <sample>_cluster_report.csv
    out_cluster.close()
    out_sample.close()
                

"""
extract_reads_output_sample
:params args
:params fastafiles = list of fasta files to record ID
:writes <all_samples>_sample_id.csv 
""" 
def extract_reads_output_sample(args, fasta_files):
    # create <all_samples>_sample_id.csv
    out_sample_name = args.dir + "/" + args.output + "_sample_id.csv"
    print("Writing sample id file to", out_sample_name,"\n")
    out_sample_name_file = open(out_sample_name,'w')
    out_sample_name_file.write("id,primer\n")
    
    # loop through each fasta file to record ID to one file
    for f in fasta_files: 
        extract_reads_output_cluster(args, f, out_sample_name_file)
    
    # close <all_samples>_sample_id.csv
    out_sample_name_file.close()


def main():
    parser = argparse.ArgumentParser(description="Generating input files for get_abundance_post_collapse.py and demux_ont_with_genome.py")
    parser.add_argument('fasta', help='\t\tDirectory of fasta files.')
    parser.add_argument('-s','--samples', nargs='+', required=False, help='\t\tSpecified samples to include. Default: all files with .fa/.fasta')
    parser.add_argument('--batch', default=False, action="store_true", help='\t\tInclude batch name into output files')
    parser.add_argument('-i','--index', required=False, default = ".fa", help='\t\tSpecific string in fasta files to capture in input directory')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output sample file documenting reads.')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of fasta files.')
    
    args = parser.parse_args()
    if args.output is None:
        args.output = "output" 
    
    if args.dir is None:
        args.dir = args.fasta 
        
    fasta_files = list_named_fasta(args)
    extract_reads_output_sample(args, fasta_files)
    
    print("All Done")
    
if __name__ == "__main__":
  main()
