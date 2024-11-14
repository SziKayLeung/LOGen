#!/usr/bin/env python


"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Replace PB.ID in gtf, gff, fasta, and classification 

Example: PB.1.2 --> ONT1.1.2 
Replacing PB with ONT<chromosome_number>

Chromosome number can be determined in gtf, gff and classification file, but not in fasta file
Therefore to replace fasta file, also needs gtf/gtf/classification input
"""

import os
from gtfparse import read_gtf
from Bio import SeqIO
import argparse


def replace_gtf(args):
    output_gtf = args.dir + str("/") + args.output_gtf
    print("Writing output to", output_gtf)

    # ID dictionary <old_ID>:<new_ID>
    id_dict = {}
    with open(args.gtf) as f, open(output_gtf, 'w') as out_file:
        for line in f:
            # Split line by tabs and strip unnecessary whitespace
            columns = [column.strip() for column in line.rstrip('\n').split('\t')]
    
            # keep the chromosome number
            replaceTxt = columns[0].replace("chr", "ONT")
    
            # Attributes column 
            replaceAttributes = columns[8].replace("PB", replaceTxt)

            # create dictionary    
            original_id = columns[8].split(";")[0].split(" ")[1].replace('"', '')
            transcript_id = replaceAttributes.split(";")[0].split(" ")[1].replace('"', '')
            id_dict[original_id] = transcript_id
            
            # keep all columnns but substitute new attributes column
            newOutput = columns[0:8] + [str(replaceAttributes)]
    
            # write output
            out_file.write('\t'.join(newOutput) + '\n') 

    return(id_dict)
    


def replace_classification(args):
    output_classfile = args.dir + str("/") + args.output_class
    print("Writing output to", output_classfile)

    with open(args.classification) as f, open(output_classfile, 'w') as out_file:
        id_dict = {}
        for i, line in enumerate(f):
            columns = [column.strip() for column in line.rstrip('\n').split('\t')]
            
            # only extract for lines after header
            if i != 0:
            	# columns 0: transcript_id
            	# columns 1: chromosome
                replaceTxt = columns[1].replace("chr", "ONT")
                
                # create dictionary
                original_id = columns[0]
                transcript_id = columns[0].replace("PB", replaceTxt)
                id_dict[original_id] = transcript_id

                new_output = [str(transcript_id)] + columns[1:] 
            else:
                new_output = columns
    
            out_file.write('\t'.join(new_output) + '\n') 

    return(id_dict)

# use the id_dict from classification or gtf/gff file
def replace_fasta(args, id_dict):
	# read fasta sequence
    fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')

    output_fasta = args.dir + str("/") + args.output_fasta
    print("Writing output to", output_fasta)

    with open(output_fasta, "w") as f:
        for seq in fasta_sequences:
            if seq.id in id_dict:
                seq.id = seq.description = id_dict[seq.id]
            else:
            	print("Note: ", seq.id, " is not in classification file, keeping original ID")
            SeqIO.write([seq], f, "fasta")


def main():
    parser = argparse.ArgumentParser(description="Subset fasta and gtf(gff) from list of ID")
    parser.add_argument('-g', '--gtf', default = None, required = False, help='gtf file to replace')
    parser.add_argument('--gff', default = None, required = False, help='gff file to replace')
    parser.add_argument('-c', '--classification', default = None, required = False, help='suffix: txt, classification txt file to replace')
    parser.add_argument('-f', '--fasta', default = None, required = False, help='suffix: fa, fasta file to replace, need either --gtf or --classification')
    parser.add_argument('-o','--output', default="renamed", required=False, help='\t\tPrefix appended for output file; Default: renamed.')
    parser.add_argument('-d','--dir', default=None, required=False, help='\t\tDirectory for output files. Default: directory of input file.')
  
    args = parser.parse_args()
    if (args.fasta is not None and args.classification is None) and (args.fasta is not None and args.gtf is None):
        print("Error: need classification or gtf file to extract the chromosome number")
        sys.exit(1)
    
    if args.gtf is not None: 
        if args.dir is None: args.dir = os.path.abspath(os.path.dirname(args.gtf))
        if args.output == "renamed":
            args.output_gtf = os.path.basename(args.gtf).replace(".gtf", "") + "_renamed.gtf"
        else:
            args.output_gtf = args.output + ".gtf"
        print("****Replacing gtf file")
        id_dict = replace_gtf(args)

    if args.classification is not None: 
        if args.dir is None: args.dir = os.path.abspath(os.path.dirname(args.classification))
        if args.output == "renamed":
            args.output_class = os.path.basename(args.classification).replace(".txt", "") + "_renamed.txt"
        else:
            args.output_class = args.output + ".txt"
        print("****Replacing classification file")
        id_dict = replace_classification(args)
    
    if args.fasta is not None: 
        extension = os.path.splitext(os.path.basename(args.fasta))[-1]
        if args.dir is None: args.dir = os.path.abspath(os.path.dirname(args.fasta))
        if args.output == "renamed":
            args.output_fasta = os.path.basename(args.fasta).replace(extension, "") + "_renamed" + extension
        else:
            args.output_fasta = args.output + extension
        print("****Replacing fasta file")
        replace_fasta(args, id_dict)

main()
