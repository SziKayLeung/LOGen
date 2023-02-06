#!/usr/bin/env python3

import pandas as pd
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description="Creating a cluster_report.csv format as input for cupcake tofu get_abundance_post_collapse.py")
    parser.add_argument("--cluster", help='\t\tTALON output annot_reads.tsv.')
    parser.add_argument('--o', help='\t\toutput dir')
    
    args = parser.parse_args()
    
    cluster = pd.read_csv(args.cluster, header = None)
    cluster.columns = ["cluster_id"]
    cluster.loc[:,"read_id"] = cluster.loc[:, 'cluster_id']
    cluster.loc[:,"cluster_id"] = 'transcript/' + cluster['cluster_id'].astype(str)
    cluster.loc[:,"read_type"] = "FL"
    
    cluster.to_csv(args.o, index = False)
    
    print("All Done")
    
if __name__ == "__main__":
  main()
    
