#!/usr/bin/env python


# packages
from gtfparse import read_gtf
import argparse
import os
import pandas as pd


def subset_reference_gtf(args, subs):
    print("Reading:", args.input)
    gtf = read_gtf(args.input) 
  
    # subset only for transcripts 
    transFeature = gtf.loc[gtf["feature"] == "transcript"]
    cols = ["gene_name","gene_id","transcript_name",'transcript_type',"transcript_id"]
    
    if args.gname:
        print("subsetting by gene name")
        outputFile = transFeature.loc[transFeature["gene_name"].isin(subs)][cols]
    
    if args.ens:
        print("subsetting by ENSEBML")
        outputFile = transFeature.loc[transFeature["gene_id"].isin(subs)][cols]
        
    outputFile.columns = cols
    
    outputName = args.dir + "/" + args.output + ".txt"
    print("Writing:", outputName)      
    outputFile.to_csv(outputName, index = False, sep = "\t")


def main():
    parser = argparse.ArgumentParser(description="Subset fasta and gtf from list of ID")
    parser.add_argument('input', help='\t\tInput file to subset: gtf, fasta or bed')
    parser.add_argument('-i', '--id_list', required=False, help='\t\tTxt file of list of Gene names/ENSEMBL ID to subset, no rownames and colnames')
    parser.add_argument('-I', '--ID_list', required=False, nargs='+', help='\t\tlist of Gene names/ENSEMBL ID  to subset - direct on command line delimited by ,')
    parser.add_argument('--gname', default=False, action="store_true", help='\t\tTo subset by gene names')
    parser.add_argument('--ens', default=False, action="store_true", help='\t\tTo subset by gene ENSEMBL ID')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output file; Default: out.')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: directory of input file.')
    
    args = parser.parse_args()
    args.basename = '.'.join(os.path.basename(args.input).split(".")[:-1])
    #print(basename)
    
    if not args.gname and not args.ens:
        print("Error: need to include --gname or --ens to define input format")
        sys.exit()
    
    if args.output is None:
        args.output = "out" 
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.input)
    
    if args.id_list is not None:
        print("Reading in", args.id_list)
        subs = pd.read_csv(args.id_list, sep = "\t", header = None)
        subs = list(subs[0].values)
    elif args.ID_list is not None:
        subs = [i.split(",") for i in args.ID_list][0]
    else:
        print("Error: need to include --id_list or --ID_list")
        
        
    subset_reference_gtf(args, subs)
    
    print("All Done")


if __name__ == "__main__":
  main()
