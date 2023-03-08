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
from gtfparse import read_gtf
import read_write as rw
import argparse
import pandas as pd
from Bio import SeqIO

"""
Subset gtf based on transcript ID 
:param retainedID: list of transcripts  
:param input_gtf: path of input gtf for subsetting
:param output name: string of name to be included in the output
:param *args[0]: optional argument, if present then generate the output file to this directory
:returns output_dict: dictionary of read name as key and sequence as value
"""
def subset_gtf(retainedID, input_gtf, name, **kwargs):
    
    print("Reading gtf:", input_gtf)
    gtf = read_gtf(input_gtf) 
    
    # variables
    basename = os.path.basename(input_gtf)
    outname = basename.replace(".gtf", "_" + name +".gtf")
    optdir = kwargs.get('dir', None)
    if optdir is not None:
        outdir = optdir
    else:
        outdir = os.path.dirname(input_gtf)
 
    subset_gtf = gtf[gtf['transcript_id'].isin(retainedID)]
    subset_gtf = subset_gtf[["seqname","source","feature","start","end","score","strand","frame","transcript_id","gene_id"]]
    
    print("Number of unique transcripts subsetted:", len(set(subset_gtf["transcript_id"].values)))
    print("Writing:", outdir + "/" + outname)  
    rw.writeGTF(subset_gtf, outdir + "/" + outname)
    
    return(subset_gtf)


"""
Subset keys of a dictionary from a list of keys 
:param all_sequence_dict: dictionary of fasta sequences - keys = id, values = sequence
:param selected_id_lst: df of id to retain and used for subsetting (1 column with no name)
:returns output_dict: dictionary of read name as key and sequence as value
"""    
def subset_seq_from_lst(all_sequence_dict, selected_id_lst):
    
    # name 1st column for later subsetting
    selected_id_lst.columns = ["ID"]
    
    # create empty list to save id of transcripts with match to human MAPT
    Transcript_id = []
    Transcript_seq = []

    # loop through each sequence in dictionary to find if match with human MAPT
    for idname, seq in all_sequence_dict.items():
        if idname in selected_id_lst[["ID"]].values:
            print("Extract fasta for:", idname)
            Transcript_id.append(idname)
            Transcript_seq.append(seq)
    
    output_dict = dict(zip(Transcript_id,Transcript_seq))
    return(output_dict)
    

"""
Subset fasta based on list of transcripts
:param retainedID: list of transcripts  
:param input_fasta: path of fasta for subsetting
:param name: string of name to be included in the output
:writes subsetted fasta
"""      
def subset_fasta(retainedID, input_fasta, name):
    fasta_sequences = SeqIO.parse(open(input_fasta),'fasta')
    
    basename = '.'.join(os.path.basename(input_fasta).split(".")[:-1])
    outputfile = basename + "_" + name +".fa"
    print("Writing to:", outputfile)
    with open(outputfile, "w") as f:
        for seq in fasta_sequences:
            if seq.id in retainedID:
                SeqIO.write([seq], f, "fasta")


def main():
    parser = argparse.ArgumentParser(description="Subset fasta and gtf from list of ID")
    parser.add_argument('input', help='\t\tInput file to subset: gtf or fasta')
    parser.add_argument('-i', '--id_list', help='\t\tTxt file of list of ID to subset, no rownames and colnames')
    parser.add_argument('--fa', default=False, action="store_true", help='\t\tTo subset fasta input file')
    parser.add_argument('--gtf', default=False, action="store_true", help='\t\tTo subset gtf input file')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output file; Default: out.')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: directory of input file.')
  
    args = parser.parse_args()
    args.basename = '.'.join(os.path.basename(args.input).split(".")[:-1])
    #print(basename)
    
    if not args.fa and not args.gtf:
        print("Error: need to include --fa or --gtf to define input format")
        sys.exit()
    
    if args.output is None:
        args.output = "out" 
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.input)
    
    print("Reading in", args.id_list)
    id_list = pd.read_csv(args.id_list, sep = "\t", header = None)
    retainedID = list(id_list[0].values)
    
    if args.fa:
        subset_fasta(retainedID, args.input, args.output)
        
    if args.gtf:
        subset_gtf(retainedID, args.input, args.output)       
   
    print("All Done")
    

if __name__ == "__main__":
    main()
