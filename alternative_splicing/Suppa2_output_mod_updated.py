#!/usr/bin/env python
# Szi Kay Leung: sl693@exeter.ac.uk
# 09/03/2020: create script to replace PB.ID in SUPPA2 output with associated gene name in SQANTI2 classification file

import pandas as pd
import os
import os.path
import sys

name=sys.argv[1]
SUPPA2_input_dir=sys.argv[2]    # input directory of all output files from SUPPA2 
SQANTI2_classification_file=sys.argv[3] # classification file (SQANTI2 filtered) 

print("Looping through dataset: " + name) 
print("Working with SUPPA2 files in: " + SUPPA2_input_dir)
print("Working with SQANTI2 classification file: " + SQANTI2_classification_file)

########### Functions 
# suppa2_input
# Aim: reads in the list of SUPPA2 output ioe files in directory and saves as list 
# Input: SUPPA2 directory containing SUPPA2 outut ioe files, name of the prefix of the ioe files 
# Output: "SUPPA2_output" iist of ioe files 
# Example: suppa2_input("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Isoseq_Paper/Human_SUPPA", "adult") 
def suppa2_input(input_SUPPA2_dir, name):
    SUPPA2_output = []
    
    for file in os.listdir(input_SUPPA2_dir):
        if file.startswith(name) and file.endswith(".ioe"):
            print(file)
            SUPPA2_output.append(os.path.join(input_SUPPA2_dir, file))
    
    return(SUPPA2_output)

# sqanti2_class_prepare
# Aim: Create dictionary of PB.id and associated_gene/associated_transcript for replacement in SUPPA2 output files
# Input: SQANTI2 classification file (filtered)
# Output: "id_dict" dictionary values
def sqanti2_class_prepare(input_SQANTI2_file):
    sqanti2 = pd.read_csv(input_SQANTI2_file, sep='\t', header = 0)
    
    # create ID for replacement 
    ID = sqanti2["associated_gene"] + '_' + sqanti2["associated_transcript"] + '_' + sqanti2["isoform"] + '_' + sqanti2["structural_category"]
    
    # Create dictionary from SQANTI2 input file of PB.ID (Isoform column from SQANTI classification file) and ID created above
    list_of_transcripts = ID.to_list()
    list_of_PB_id = sqanti2['isoform'].to_list()
    id_dict = dict(zip(list_of_PB_id,list_of_transcripts))
    
    return(id_dict)

# replace_words
# Aims: With the input value, finds the key and replaces with value from dictionary and saves to new file
# https://overlaid.net/2016/02/08/replace-words-in-files-or-strings-using-python/
def replace_words(input_file, device_values):
    file = open(input_file, 'r')
    filetemp = file.read()
    file.close 
    
    for key, val in device_values.items():
        filetemp = filetemp.replace(key, val)
        
    return filetemp


# loop_replace_words_suppa2
def loop_replace_words_suppa2(SUPPA2_output, id_dict):
    # Parse through files in SUPPA2_output from suppa2_input function
    for file in SUPPA2_output: 
        output_name = os.path.basename(file)
        print("Processing with:" + file)
        
        # Apply replace_words function with id_dict from sqanti2_class_prepare function 
        result = replace_words(file, id_dict)
        
        # Create new file and save result 
        output_file = str(output_name) + ".txt"
        fout = open(output_file, 'w')
        fout.write(result)
        fout.close

def suppa2_mod_replacement(SUPPA2_input_dir, SQANTI2_input_file):
    dat_SUPPA2_output = suppa2_input(SUPPA2_input_dir, name)
    dat_id_dict = sqanti2_class_prepare(SQANTI2_input_file)
    os.chdir(SUPPA2_input_dir)
    loop_replace_words_suppa2(dat_SUPPA2_output,dat_id_dict)
    
########### Apply Function 
suppa2_mod_replacement(SUPPA2_input_dir, SQANTI2_classification_file)
