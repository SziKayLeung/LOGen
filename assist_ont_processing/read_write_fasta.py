#!/usr/bin/env python


"""
Szi Kay Leung
S.K.Leung@exeter.ac.uk

Read and write fasta file
"""


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


"""
read input.fasta and save as dictionary(pbid:sequence)
params::
    inputfasta = input fasta file with header
return::
    output_dict = dictionary of read name as key and sequence as value
"""
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



"""
prepare_output: write dictionary to folder
params::
    subset_dict = subsetted fasta (as dictionary: keys - readname, values - sequence)
    outputdir = output folder
    outputname = name of file (.fasta will be appended)
"""
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
