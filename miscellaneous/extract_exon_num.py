#!/usr/bin/env python

# Aim: parse through a gtf to output summary gene and transcript statistics
# rough running time for gtf 1.4Kb = 11minutes
# output: csv file <ENSEMBLID>;<GENENAME>;<GENETYPE>;<TRANSCRIPT>;<NUMEXON>

# packages
import os
import argparse
import gffutils
import pandas as pd
import csv


# extract_exon_num 
def extract_exon_num(args):
    # create gffutils database
    print("Create database")
    db = gffutils.create_db(args.gtf, dbfn="input.db", force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)
    
    with open(args.dir + args.output, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        print("Writing output to:", args.dir + args.output)
        csv_writer.writerow(['Gene', 'GeneName', 'GeneType', 'Transcript', 'NumExon'])

        for transcript in db.features_of_type('transcript'):
            print("Extracting:", transcript.id)

            # extract gene features for each transcript looped
            gene = transcript.attributes['gene_id'][0]
            genename = transcript.attributes['gene_name'][0]
            gtype = transcript.attributes['gene_type'][0]

            # for each transcript, extract the number of exons by looping 
            exon_count = sum(1 for _ in db.children(transcript, featuretype='exon', order_by='start'))

            # output tuple
            #output.append((gene, genename, gtype, transcript.id, exon_count))
            csv_writer.writerow([gene, genename, gtype, transcript.id, exon_count])
            
    db.conn.close()

def main():
    parser = argparse.ArgumentParser(description="Extract the number of exons for each transcript in gtf")
    parser.add_argument('gtf', help='\t\tInput gtf')    
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for output file; Default: XX.numExon.csv')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: directory of input file.')
    
    args = parser.parse_args()
    
    if args.dir is None:
        args.dir = os.path.dirname(args.gtf) or os.getcwd()
    
    if args.output is None: 
        args.basename = '.'.join(os.path.basename(args.gtf).split(".")[:-1])
        args.output = "/" + args.basename + ".numExon.csv"
    else:
        args.output = "/" + args.output + ".csv"
    
    extract_exon_num(args)
    print("All done")

if __name__ == "__main__":
    main()
