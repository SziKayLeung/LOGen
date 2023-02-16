#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Date: April 2022
Aim: fill in itemRgb column in bed file according to abundance and protein-coding potential
# colour reference: https://htmlcolorcodes.com/colors/shades-of-green/

---- Input:
==> Abundance: 
# if abundance file contains multiple datasets
isoform,dataset1_sample1,dataset2_sample2...,dataset1_sum_FL,dataset2_sum_FL,dataset
PB.1.1 ,1               ,2 ........         ,20             ,40             ,Both
PB.2.1 ,0               ,2 ........         ,0              ,20             ,dataset1
PB.3.1 ,2               ,0 ........         ,20             ,0              ,dataset2

# OR
isoform,MergedFL
PB.1.1,60
PB.1.2,20

"""


# packages
import pandas as pd
import numpy as np
import argparse
import os

# turn off warnings in this script (relating to loc)
pd.options.mode.chained_assignment = None  # default='warn'


"""
read in files for downstream
abundance as dataframe
coding potential cut-off dependent on species (https://cpat.readthedocs.io/en/latest/)
:return list:
        counts = read in abundance df
        bed = read in bed
        coding = dictionary; {isoform, Coding_Status <Coding, Non_Coding, No ORF>}
                                  {isoform, Coding_Status <None>} if not using cpat input files
""" 
def input_files(args):
    
    # read in abundance file
    print("Reading abundance file")
    abundance = pd.read_csv(args.a)
    
    # read in bed file 
    print("Reading bed file")
    bed = pd.read_csv(args.bed, sep = "\t", header = None)
    bed.columns = ["chr","start","end","transcript","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
    
    if args.cpat is not None:
        print("Reading CPAT input files")
        cpat = pd.read_csv(args.cpat, sep = "\t")
        cpat_noORF = pd.read_csv(args.noORF, sep = "\t", header = None)
        
        # default cut off thresholds
        if args.species == "human": 
            threshold = 0.364 
        else:
            threshold = 0.44
        print("Using", threshold,"as cut-off for coding potential")
        
        # determine coding status based on threshold
        cpat["Coding_Status"] = np.where(cpat["Coding_prob"] >= threshold, "Coding", "Non_Coding")
        cpat_noORF["Coding_Status"] = "No_ORF"
        
        coding_dict = {**dict(zip(cpat["seq_ID"],cpat["Coding_Status"])),**dict(zip(cpat_noORF[0],cpat_noORF["Coding_Status"]))}
    else:
        coding_dict = {keys:"None" for keys in abundance["isoform"]}
    
    # wrap into one list for output    
    inputs = {
        'counts' : abundance,
        'bed' : bed,
        'coding' : coding_dict,
    }
    
    return inputs

    
"""
colour_palette
determine the colour based on abundance and protein coding potential
:params abundance = int: full-length read count
:params status = str: protein-coding potential <Non_Coding, Coding, NA>
:return colour = list of RGB 
""" 
def colour_palette(abundance, status):    
    if status == "Non_Coding":
        FL25 = [250,160,160] # pastel red
        FL50 = [227,83,53] # poppy
        FL100 = [70,74,68] # brick red
        FL250 = [255,0,0] # red
        FLmore = [238,75,43] # bright red
    elif status == "Coding":
        FL25=[193,225,193] # pastel green
        FL50 = [147,197,114] # pistachio
        FL100 = [80,200,120] # emerald green
        FL250 = [34,139,34] # forest green
        FLmore = [0,128,0] # green
    else: 
        FL25=[211,211,211] # pastel grey
        FL50 = [192,192,192] # silver
        FL100 = [169,169,169] # dark grey
        FL250 = [113,121,126] # steel grey
        FLmore = [0,0,0] # black
    
    if abundance <= 25:
        return(FL25)
    elif abundance > 25 and abundance <= 50:
        return(FL50)
    elif abundance > 50 and abundance <= 100:
        return(FL100)
    elif abundance > 100 and abundance <= 250:
        return(FL250)
    else:
        return(FLmore) 

    return(col)


"""
colour the bed file per row (isoform) based on the abundance and coding potential
apply the colour_palette()
the output bed file only keeps the isoforms in the abundance dictionary
:params inputs = list output from input_files()
:params abundance_dict = dictionary of {isoform: read_count} 
:params output_name = str of name for file 
:writes coloured bed file
    1. _uniqueDataset1counts_coloured.bed12
    1. _bothDataset1counts_coloured.bed12
""" 
def colour_uscs_tracks(inputs, abundance_dict, output_name):
    
    # subset the bed file by only keeping the trascripts that are in the abundance dictionary
    # enables multiple bed files generated based on the input abundance dictionary
    outputBed = inputs['bed'].loc[inputs['bed']["transcript"].isin(abundance_dict.keys())]
    
    # loop through each row of subsetted bed file and determine the colour
    col_output = []
    for index, row in outputBed.iterrows():
        transcript = row[3] 
        a = int(abundance_dict[transcript])
        c = inputs['coding'][transcript] 
        if c == "Coding":
            col_output.append(colour_palette(a, "Coding"))
        elif c == "Non_Coding":
            col_output.append(colour_palette(a, "Non_Coding"))
        else:
            col_output.append(colour_palette(a, "Non_ORF"))
    
    # replace the colour output list to itemRgb colour 
    # remove strings [,] (output of colour_palette())
    outputBed['itemRgb'] = col_output
    outputBed['itemRgb'] = [str(i).replace("[", "") for i in outputBed['itemRgb']]
    outputBed['itemRgb'] = [str(i).replace("]", "") for i in outputBed['itemRgb']]
    outputBed['itemRgb'] = [str(i).replace(" ", "") for i in outputBed['itemRgb']]
    
    basename = os.path.basename(output_name)
    print("Writing to:", basename)
    outputBed.to_csv(output_name, sep = "\t", header = None, index = None)


"""
generate multiple coloured bed files per dataset (if abundance file contains 'dataset' column)
:params inputs = list output from input_files()
:params datasetName = string of different datasetname in "dataset" column in abundance file
:writes 2 coloured bed file per dataset
""" 
def col_by_dataset(args, inputs, datasetName):
    
    '''
    Generates two bed files per dataset (i.e ONT, PacBio)
    if datasetName = "ONT"
    1) common isoforms that are detected in both datasets (i.e. ONT and PacBio), but using the ONT FL read counts
    2) unique isoforms that are only detected in ONT dataset 
    
    Pre-requisites for abundance file:
    - dataset column, containing the datasetName (i.e. Both, ONT, PacBio)
    - sum FL columns for each dataset, containing the datasetName (i.e. ONT_sum_FL, PacBio_sum_FL)
    '''
    
    print("Processing", datasetName)
    
    ## --- Both common isoforms
    # obtain common isoforms that are present in both datasets
    # create an abundance dictionary, with read counts of these common isoforms from dataset of interest
    # subset bed file using abundance dictionary (isoforms of interest as key) and colour
    bothdat = inputs['counts'].loc[inputs['counts']["dataset"] == "Both"]
    both_a_dict = dict(zip(bothdat["isoform"],bothdat[datasetName + "_sum_FL"]))
    bothoutname = args.dir + args.o + "_both" + str(datasetName) + "counts_coloured.bed12"
    colour_uscs_tracks(inputs, both_a_dict, bothoutname)
    
    ## --- unique isoforms
    # obtain unique isoforms to dataset
    # obtain the read counts of these unique isoforms
    unidat = inputs['counts'].loc[inputs['counts']["dataset"] == datasetName]  
    unidat_a_dict = dict(zip(unidat["isoform"],unidat[datasetName + "_sum_FL"]))
    unioutname = args.dir  + args.o + "_unique" + str(datasetName) + "counts_coloured.bed12"
    colour_uscs_tracks(inputs, unidat_a_dict, unioutname)


"""
generate one coloured bed file (if abundance file contains 'dataset' column)
read counts are concatenated across datasets
:params inputs = list output from input_files()
:writes one coloured bed file
    _concat_counts_coloured.inputs.bed12
""" 
def col_by_all(args, inputs):
    
    '''
    Generates one bed file with counts concatenated
    common isoforms that are detected in both datasets (i.e. ONT and PacBio), sum reads counts from "sum" columns
    unique isoforms that are only detected in 1 dataset, include only counts from that dataset

    same pre-requisites for abundance file col_by_dataset()
    '''
    
    # filer on columns with "sum" in name (i.e ONT_sum_FL, PacBio_sum_FL)
    # for unique isoforms, 0 reads in one column
    inputs['counts']["all_sum_FL"] = inputs['counts'].filter(regex='sum').sum(axis=1)
    
    # create abundance dictionary and colour
    abundance_dict = dict(zip(inputs['counts']["isoform"],inputs['counts']["all_sum_FL"]))
    concatname = args.dir + args.o + "_concat_counts_coloured.inputs.bed12"
    colour_uscs_tracks(inputs, abundance_dict, concatname)
    
    
def main():
    parser = argparse.ArgumentParser(description="Colour Transcripts from ONT and Iso-Seq Targeted Transcriptome by abundance and coding potential")
    parser.add_argument('--bed', "--input_bed", help='\t\tSorted bed12 file from the final merged dataset')
    parser.add_argument('--cpat', "--cpat_ORF", help='\t\tCPAT output file of best ORFs')
    parser.add_argument('--noORF',"--CPAT_noORF", help='\t\tCPAT output file of transcripts with no ORFs')
    parser.add_argument('--a',"--abundance", help='\t\tAbundance file of FL reads from final merged dataset')
    parser.add_argument('--o',"--output", help='\t\tOutput coloured bed12 file from the final merged dataset')
    parser.add_argument('-s','--species', default = None, choices=['human', 'mouse'], required=False, help='\t\tDirectory for output files. Default: Directory of input bed.')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of input bed.')
    
    args = parser.parse_args() 
    
    #---- QC arguments
    
    if args.dir is None:
        args.dir = os.path.dirname(args.bed) + "/"
    print("Writing output files to:", args.dir)
    
    # required input if input coding potential file   
    if args.cpat is not None:
        if args.noORF is None:
            print("Error: CPAT output file of transcripts with no ORFs required")
        if args.species is None:
            print("Error: --species required: mouse, human")

    #---- run  
    
    inputs = input_files(args)
    
    # if dataset is in the column, then run col_by_dataset() and col_by_all()
    if 'dataset' in inputs['counts'].columns:
        
        # generate one merged bed file with concatenated counts
        col_by_all(args, inputs)
        
        # generate 2 bed files for each dataset 
        for datasetName in set(inputs['counts']["dataset"].values):
            if datasetName != "Both":
                col_by_dataset(args, inputs, datasetName)
    
    else:
        abundance_dict = dict(zip(inputs['counts']["isoform"],inputs['counts']["MergedFL"]))
        colour_uscs_tracks(inputs, abundance_dict, args.dir + args.o + "_coloured.bed12")
        
if __name__ == "__main__":
    main()
