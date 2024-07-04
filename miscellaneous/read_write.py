#!/usr/bin/env python

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Read and write fasta/gtf file

Functions:
    readFasta(inputfasta)
    writeFasta(subset_dict, outputdir, outputname)
    writeGTF(inGTF,file_path)
"""

import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


"""
Read input.fasta and save as dictionary(pbid:sequence)
:param inputfasta: input fasta file with header 
:returns output_dict: dictionary of read name as key and sequence as value
"""
def readFasta(inputfasta):
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
Write fasta directory (key = isoform id; value = fasta sequence) to a fasta file 
:param subset_dict: subsetted fasta (as dictionary: keys - readname, values - sequence) 
:param outputdir: path/to/the/file.fasta
:param outputname: name of file (.fasta will be appended)
:returns nothing
"""
def writeFasta(subset_dict, outputdir, outputname):
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


"""
Write a GTF dataframe into a file
:param inGTF: GTF dataframe to be written. It should either have 9 columns with the last one being the "attributes" section or more than 9 columns where all columns after the 8th will be colapsed into one.
:param file_path: path/to/the/file.gtf
:returns nothing
https://github.com/mpg-age-bioinformatics/AGEpy/blob/master/AGEpy/gtf.py
"""
def writeGTF(inGTF, file_path):
    cols = inGTF.columns.tolist()
    if len(cols) == 9:
        print("yes")
        if 'attribute' in cols:
            df = inGTF.copy()
    else:
        df = inGTF[cols[:8]].copy()
        df.loc[:, 'attribute'] = ""
        for c in cols[8:]:
            if c == cols[len(cols) - 1]:
                df.loc[:, 'attribute'] += c + ' "' + inGTF[c].astype(str) + '";'
            else:
                df.loc[:, 'attribute'] += c + ' "' + inGTF[c].astype(str) + '"; '

    df.to_csv(file_path, sep="\t", header=None, index=None, quoting=csv.QUOTE_NONE)


"""
Write list to output txt file
:param outtxt: path/to/output/file.txt
:param lst: list of elements to write out line by line
:returns nothing
"""    
    
def write_lst(outtxt, lst):
    
    textfile = open(outtxt, "w")
    print("Writing output to:", outtxt)
    for element in lst:
        textfile.write(element + "\n")
    textfile.close()
