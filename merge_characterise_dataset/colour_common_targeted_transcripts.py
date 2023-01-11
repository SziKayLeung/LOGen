#!/usr/bin/env python3
# Szi Kay Leung (sl693@exeter.ac.uk)
# colour reference: https://htmlcolorcodes.com/colors/shades-of-green/

import pandas as pd
import numpy as np
import argparse

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

def colour_uscs_tracks(args):
    bed = pd.read_csv(args.bed, sep = "\t", header = None)
    
    print("**** Reading abundance files.")
    abundance = pd.read_csv(args.a, sep = ",")
    abundance_dict = dict(zip(abundance["id"],abundance["MergedFL"]))
    
    if args.cpat is not None:
        print("**** Reading CPAT input files.")
        cpat = pd.read_csv(args.cpat, sep = "\t")
        cpat_noORF = pd.read_csv(args.noORF, sep = "\t", header = None)
        
        cpat["Coding_Status"] = np.where(cpat["Coding_prob"] >= 0.44, "Coding", "Non_Coding")
        cpat_noORF["Coding_Status"] = "No_ORF"
        
        coding_dict = {**dict(zip(cpat["seq_ID"],cpat["Coding_Status"])),**dict(zip(cpat_noORF[0],cpat_noORF["Coding_Status"]))}
    else:
        coding_dict = {keys:"None" for keys in abundance["id"]}
    
    col_output = []
    for index, row in bed.iterrows():
        transcript = row[3] 
        a = int(abundance_dict[transcript])
        c = coding_dict[transcript] 
        if c == "Coding":
            col_output.append(colour_palette(a, "Coding"))
        elif c == "Non_Coding":
            col_output.append(colour_palette(a, "Non_Coding"))
        else:
            col_output.append(colour_palette(a, "Non_ORF"))
    
    bed[8] = col_output
    bed[8] = [str(i).replace("[", "") for i in bed[8]]
    bed[8] = [str(i).replace("]", "") for i in bed[8]]
    bed[8] = [str(i).replace(" ", "") for i in bed[8]]
    
    bed.to_csv(args.o, sep = "\t", header = None, index = None)

def main():
    parser = argparse.ArgumentParser(description="Colour Transcripts from ONT and Iso-Seq Targeted Transcriptome by abundance and coding potential")
    parser.add_argument('--bed', "--input_bed", help='\t\tSorted bed12 file from the final merged dataset')
    parser.add_argument('--cpat', "--cpat_ORF", help='\t\tCPAT output file of best ORFs')
    parser.add_argument('--noORF',"--CPAT_noORF", help='\t\tCPAT output file of transcripts with no ORFs')
    parser.add_argument('--a',"--abundance", help='\t\tAbundance file of FL reads from final merged dataset')
    parser.add_argument('--o',"--output", help='\t\tOutput coloured bed12 file from the final merged dataset')
    
    args = parser.parse_args() 
    colour_uscs_tracks(args)

if __name__ == "__main__":
    main()
