#!/usr/bin/env python3

import pandas as pd
import csv
import argparse
import glob
import os


def read_id(filename):
  print(filename)
  df = pd.read_csv(filename, sep = "\t", header = None) 
  df.columns = ["id"]
  sample = os.path.basename(filename).replace("_clean_id.txt","")
  df["primer"] = sample

  return(df)


def main():
  parser = argparse.ArgumentParser(description="Creating a cluster_report.csv format as input for cupcake tofu get_abundance_post_collapse.py")
  parser.add_argument("--dir", help='\t\tTALON output annot_reads.tsv.')
  parser.add_argument('--o', help='\t\toutput dir')

  args = parser.parse_args()
  
  files = [read_id(i) for i in glob.glob(args.dir + '/**/**_clean_id.txt', recursive=True)]
  filesdf = pd.concat(files)
  
  filesdf.to_csv(args.o, index = False)
  
  print("All Done")

if __name__ == "__main__":
  main()

