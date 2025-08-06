#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def convert_nucleotide_to_protein(args):
  dir = os.path.dirname(args.fasta) + "/"
  name = os.path.basename(args.fasta).split(".")[0]
  print(dir + name + "protein.fasta")
  with open(args.fasta) as infile, open(dir + name + "protein.fasta", "w") as outfile:
      for record in SeqIO.parse(infile, "fasta"):
          protein_seq = Seq(record.seq).translate()
          outfile.write(f">{record.id}_translated\n{protein_seq}\n")


def main():
    parser = argparse.ArgumentParser(description="Generating input files for get_abundance_post_collapse.py and demux_ont_with_genome.py")
    parser.add_argument('fasta', help='\t\tinput path fasta.')
    
    args = parser.parse_args()
    convert_nucleotide_to_protein(args)
    
if __name__ == "__main__":
  main()
