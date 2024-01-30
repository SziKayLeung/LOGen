#!/usr/bin/env python3
# Szi Kay Leung (sl693@exeter.ac.uk)
# 21/01/2022: Extract fasta for the best ORF determined using CPAT

# load packages
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import argparse

# read input.fasta and save as dictionary(pbid:sequence)
def read_inputfasta(inputfasta):
    f = open(inputfasta, "r")
    # loop through the input file and save the name and sequence using SeqIO.parse
    seq = []
    input_genes = []
    for seq_record in SeqIO.parse(f, "fasta"):
        input_genes.append(seq_record.name)
        seq.append(str(seq_record.seq))
    
    # Key: PB_ID, sequence
    output_dict = dict(zip(input_genes,seq))
    return(output_dict)

def identify_best_ORF_fasta(all_sequence_dict, input_best_ORF):
    # create empty list to save id of transcripts with match to human MAPT
    Transcript_id = []
    Transcript_seq = []

    # loop through each sequence in dictionary to find if match with human MAPT
    print("Extracting ORF from best sequence")
    for orf_id, seq in all_sequence_dict.items():
        if orf_id in input_best_ORF[["ID"]].values:
            print("Extract fasta for:", orf_id)
            Transcript_id.append(orf_id)
            Transcript_seq.append(seq)
    
    Transcript = dict(zip(Transcript_id,Transcript_seq))
    return(Transcript)

# output the result of subset fasta (now as dictionary) as a format for output
def prepare_output(subset_dict, outputdir, outputname):
    # need to save the dictionary back to list of id and sequence
    subset_id = list(subset_dict.keys())
    subset_seq = list(subset_dict.values())

    # loop through the dictionary and save each record
    link_record = []
    for count in range(len(subset_dict)):
        record = subset_id[count],subset_seq[count]
        link_record.append(record)

    # final output for SeqIO
    record_list = []
    for (pb_id, sequence) in link_record:
        record = SeqRecord(Seq(sequence),id=pb_id)
        record_list.append(record)
    
    # output IDs only to text file 
    with open(outputdir + "/" + outputname + "_Ids.txt", "w") as f:
        for key in subset_dict:
            print(key, file=f)
    
    # output fasta sequence
    SeqIO.write(record_list, outputdir + "/" + outputname + ".fasta", "fasta")


def main():
    parser = argparse.ArgumentParser(description="Extracting the fasta sequence for the best ORF predicted from CPAT")
    parser.add_argument('--fa', "--input_fasta", help='\t\tFasta sequence generated from CPAT.')
    parser.add_argument('--orf', "--best_orf", help='\t\t File generated from CPAT tabulating the best ORF only')
    parser.add_argument('--o_name',"--output_name", help='\t\t Output name.')
    parser.add_argument('--o_dir', "--output_dir", help='\t\tOutput path and name for list of ONT retained transcript IDs')

    args = parser.parse_args()
    print("Reading in ORF")
    ORF = pd.read_csv(args.orf, sep = "\t")
    all_sequence_dict = read_inputfasta(args.fa)
    
    subset_dict = identify_best_ORF_fasta(all_sequence_dict, ORF)
    prepare_output(subset_dict, args.o_dir, args.o_name)    
    print("All Done")
    
if __name__ == "__main__":
    main()

