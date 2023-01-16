#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Subset sqanti classification file by target genes using associated_gene column
Functions:
    read_and_subset_target_genes(TargetGene, sqanti_input_file)
"""

# packages
import pandas as pd
from gtfparse import read_gtf
import csv
import sys
import argparse
import os.path


"""
read and subset sqanti classification files
:param TargetGene: txt file of list of target genes for filtering; column 1: ENSEMBL ID; column 2: genename
:param sqanti_input_file: path/to/sqanti/classification/file.txt
:returns sq: df of sqanti file of target genes
:write <name>_targetgenes.txt
"""
def read_and_subset_target_genes(TargetGene, sqanti_input_file):
  
  # generate path for output file
  basename = os.path.basename(sqanti_input_file)
  outname = basename.replace(".txt","_targetgenes.txt")
  outdir = os.path.dirname(sqanti_input_file)
  
  print("Reading:", sqanti_input_file)
  df = pd.read_csv(sqanti_input_file, sep = "\t",dtype={"bite": "string","polyA_motif":"string"})
  
  print("Reading:", TargetGene)
  TargetGene_file = pd.read_csv(TargetGene, sep = "\t", header = None)
  TargetGene = list(TargetGene_file[[1]].values.ravel())
  TargetGene_Ensembl = list(TargetGene_file[[0]].values.ravel())
  
  print("Filtering with genes:")
  print(', '.join(TargetGene))
  # capitalise to ensure correct subset
  TargetGene = [x.upper() for x in TargetGene]
  target_df = df[df["associated_gene"].str.upper().isin(TargetGene)]
    
  # capture any fusion transcripts where the target gene is part of the associated_gene
  fusion = df[df["structural_category"] == "fusion"]
  target_index = []
  for index, row in fusion.iterrows():
    for gene in TargetGene:
      if gene in row["associated_gene"].upper():
        #print("Fusion transcript with target gene:", row["associated_gene"].upper())
        target_index.append(index)
  
  # rbind
  target_fusion = fusion.loc[target_index]
  df = pd.concat([target_df, target_fusion])
  
  print("Writing:", outdir + "/" + outname)
  df.to_csv(outdir + "/" + outname, sep = "\t")
  
  return(df)
  
  
def main():
    parser = argparse.ArgumentParser(description="Subset sqanti classification file by target genes")
    parser.add_argument("--target_gene", help='\t\t Txt file of list of target genes for filtering; column 1: ENSEMBL ID; column 2: genename')
    parser.add_argument('--class_file', help='\t\t SQANTI input classification file')
    
    args = parser.parse_args()
    read_and_subset_target_genes(args.target_gene, args.class_file)
 
    print("All Done")
    
if __name__ == "__main__":
    main()

