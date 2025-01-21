#!/usr/bin/env python


"""
Szi Kay Leung
S.K.Leung@exeter.ac.uk

Read and write fasta file
"""


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os.path



"""
read input.fasta and save as dictionary(pbid:sequence)
params::
    inputfasta = input fasta file with header
return::
    output_dict = dictionary of read name as key and sequence as value
"""
def read_inputfasta(inputfasta):
    if os.path.basename(inputfasta).endswith(('.fa', '.fasta')):
        print("Reading in fasta file")
        formatFasta = "fasta"
    else:
        print("Reading in fastq file")
        formatFasta = "fastq"
    f = open(inputfasta, "r")
    # loop through the input file and save the name and sequence using SeqIO.parse
    seq = []
    input_genes = []
    for seq_record in SeqIO.parse(f, formatFasta):
        input_genes.append(seq_record.name)
        seq.append(str(seq_record.seq))

    # Key: PB_ID, sequence
    output_dict = dict(zip(input_genes,seq))
    return(output_dict)



"""
prepare_output: write dictionary to folder
params::
    input_fasta = input fasta used to read_write
    subset_dict = subsetted fasta (as dictionary: keys - readname, values - sequence)
    outputdir = output folder
    outputname = name of file (.fasta will be appended)
"""
def prepare_output(inputfasta, subset_dict, outputdir, outputname):
        
    # need to save the dictionary back to list of id and sequence
    subset_id = list(subset_dict.keys())
    subset_seq = list(subset_dict.values())

    if os.path.basename(inputfasta).endswith(('.fa', '.fasta')):
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
    
        # output fasta sequence
        SeqIO.write(record_list, outputdir + "/" + outputname + ".fasta", "fasta")
         
    else:
    
        with open(outputdir + "/" + outputname + ".fastq", "w") as output_handle:
            pass  # This clears the file by opening it in write mode and immediately closing it


        with open(outputdir + "/" + outputname + ".fastq", "a") as output_handle:  
            for seq_record in SeqIO.parse(inputfasta, "fastq"):
                if seq_record.name in subset_id:
                    SeqIO.write(seq_record, output_handle, "fastq")
                
    # output IDs only to text file
    with open(outputdir + "/" + outputname + "_Ids.txt", "w") as f:
        for key in subset_dict:
            print(key, file=f)