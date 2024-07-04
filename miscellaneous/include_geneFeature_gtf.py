#!/usr/bin/env python

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Include gene feature in gtf if missing
"""

##------------- Load modules -----------

from gtfparse import read_gtf
import pandas as pd
import csv
import argparse
import os
import sys
import logging 
import warnings
import read_write as rw

warnings.simplefilter(action='ignore', category=FutureWarning) #gtfparse future warnings

logger = logging.getLogger()


##------------- Functions -----------

"""
prepare_gtf: reads gtf and includes gene name from classification file if included in paramaters
:params args
    args.classfile = <isoform> <associated_gene>; ensure isoform matches with gtf transcript_id
:returns gtf
""" 
def prepare_gtf(args):
    print("**** Reading gtf", file=sys.stdout)
    logger.disabled = True
    gtf = read_gtf(args.input_gtf)
    logger.disabled = False
    
    if args.classfile is not None:
        print("**** Reading classfile", file=sys.stdout)
        classfile = pd.read_csv(args.classfile, sep = "\t")
        
        print("**** Including genename", file=sys.stdout)
        # create a dictionary of isoform and gene name
        nameReplaceDict = dict(zip(classfile["isoform"],classfile["isoform"] + ";" + classfile["associated_gene"]))
        # replace using dictionary
        gtf['gene_id'] = gtf['transcript_id'].astype(str).map(nameReplaceDict)
        # split column and expand to include gene_name
        gtf[['gene_id','gene_name']] = gtf['gene_id'].str.split(';',expand=True)
        gtf = gtf[["seqname","source","feature","start","end","score","strand","frame","gene_id","gene_name","transcript_id"]]
    
    return(gtf)


"""
extract_gene_features: takes gtf and summarises gene feature information based on the transcript level
:params gtf from prepare_gtf
:returns dataframe of the gene-level information
""" 
def extract_gene_features(gtf):
    # extract transcript feature only - no need for exons
    tgtf = gtf.loc[gtf["feature"] == "transcript"]
    
    # loop through transcript and extract information
    genes = {}    
    for index, row in tgtf.iterrows():
        seqname = row[0]
        source = row[1]
        feature_type = row[2]
        start = int(row[3])
        end = int(row[4])
        strand = row[6]
        gene_id = row[8]
        gene_name = row[9]

        # Extract gene_id from the attributes column
        # note default source = PacBio
        if gene_id not in genes:
            genes[gene_id] = {'seqname': seqname, 
                              'source': source, 
                              'feature': 'gene', 
                              'start': start, 
                              'end': end, 
                              'strand': strand,
                              'frame': 0,
                              'gene_id': gene_id,
                              'gene_name': gene_name,
                              'transcript_id': 'NaN'}
        # ensure to encapsulate the whole range of all the transcripts associated with gene
        else:
            if start < genes[gene_id]['start']:
                genes[gene_id]['start'] = start
            if end > genes[gene_id]['end']:
                genes[gene_id]['end'] = end
    
    # transpose to dataframe and re-index
    genesDF = pd.DataFrame(genes).T
    genesDF = genesDF.reset_index(drop=True)
    
    return(genesDF)


"""
merge_finalise: merge gene-level information with initial gtf and write merged gtf
:args, gtf and genesDF from output of other two functions
:writes merged.gtf
""" 
def merge_finalise(args, gtf, genesDF):
    # merge and sort
    merge = pd.concat([gtf, genesDF])
    merge = merge.sort_values(by=['seqname','start',"transcript_id"])
    
    # write output
    logger.disabled = True
    rw.writeGTF(merge,args.output_dir + "/" + args.output_name)
    logger.disabled = False
    
    
def main():
    parser = argparse.ArgumentParser(description="Include gene feature in gtf if missing")
    parser.add_argument("-g","--input_gtf", help='\t\tInput gtf')
    parser.add_argument("-c","--classfile", help='\t\tSQANTI classification file', required=False,default=None)
    parser.add_argument("-o","--output_dir", help='\t\tOutput path, default is current working directory', required=False,default=None)
    parser.add_argument("-n","--output_name", help='\t\tOutput name, default is basename of input gtf with geneIncluded.gtf', required=False, default=None)
    
    args = parser.parse_args()
    if args.output_name is None:
        args.output_name = os.path.splitext(os.path.basename(args.input_gtf))[0] + "_geneIncluded.gtf"
    else:
        args.output_name = args.output_name + "_geneIncluded.gtf"
    
    if args.output_dir is None:
        args.output_dir =  os.getcwd()
        
    gtf = prepare_gtf(args)
    genesDF = extract_gene_features(gtf)
    merge_finalise(args, gtf, genesDF)
       
if __name__ == "__main__":
    main()
