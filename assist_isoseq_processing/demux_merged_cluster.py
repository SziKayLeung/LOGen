#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: demultiplex merged cluster_report.csv using sample flnc read id
Premise:
Individual iso-Seq datasets merged into one dataset after Iso-Seq3 refine 
Merged dataset is then clustered together (Iso-Seq3 cluster), with ccs reads across samples merged into one transcript
Script developed to demultiplex cluster_report.csv generated from Iso-Seq3 cluster, 
and retrieve origin of ccs reads clustered per sample

---- Input:
=> S1.flnc.report.csv
id,strand,fivelen,threelen,polyAlen,insertlen,primer
mxxx/ccs,..
myyy/ccs,..

=> S2.flnc.report.csv
id,strand,fivelen,threelen,polyAlen,insertlen,primer
maaa/ccs,..
mbbb/ccs,..

=> merged.clustered.cluster_report.csv
cluster_id,read_id,read_type
transcript/0,mxxx/ccs,FL
transcript/0,myyy/ccs,FL
transcript/0,maaa/ccs,FL
transcript/1,mbbb/ccs,FL

---- Output:
=> merged.flnc.report.csv
id,sample
mxxx/ccs,S1
myyy/ccs,S1
maaa/ccs,S2
mbbb/ccs,S2

==> out_clustered.demuxed.cluster_report.csv
cluster_id,id,sample
transcript/0,mxxx/ccs,S1
transcript/0,myyy/ccs,S1
transcript/0,maaa/ccs,S2
transcript/1,mbbb/ccs,S2

"""

# packages
import pandas as pd
import os
import csv
import argparse
import glob


"""
merge_flnc_reads
sample column = <sample>.flnc.report.csv
:params args
:writes demux_sample_id.csv
:return out_sample_name.df = demux_sample_id.csv as dataframe
""" 
def merge_flnc_reads(args):
    '''
    1. Read in all *flnc.report.csv in refine directory
    2. Merge all IDs into one file, using prefix of flnc.repor.csv as identifier 
    '''    
    # set up output file
    out_sample_name = args.dir + "/" + args.output + "_sample_id.csv"
    print("Writing sample id file to", out_sample_name,"\n")
    out_sample_name_file = open(out_sample_name,'w')
    out_sample_name_file.write("id,primer\n")
    
    # loop through each flnc.report.csv
    for file in glob.glob(args.refine_dir + "/*flnc.report.csv"):
        # obtain the sample name using the first part split by "."
        sample = os.path.basename(file).split(".")[0]
        with open(file, 'r') as f:
            print("Recording reads from:", file,)
            # skip header line
            next(f)
            # record each read from the file and write to output
            for line in f:
                ccsread = line.split(",")[0]
                out_sample_name_file.write(str(ccsread) + "," + sample + "\n")
    out_sample_name_file.close()
    
    # read in file out for downstream
    out_sample_name_df = pd.read_csv(out_sample_name)
    
    return(out_sample_name_df)
    

"""
demux_isoseq_cluster
:params args
:writes out_clustered.demuxed.cluster_report.csv
""" 
def demux_isoseq_cluster(args):
    '''
    1. Read in flnc reads across all samples into one file for identifying origin of ccs reads
    1. Read in cluster report from merging all samples, and rename column to ease merging
    2. Merge in flnc_reads and cluster_report by "ID" column
    '''
    
    # read in flnc reads and cluster report
    flnc_reads = merge_flnc_reads(args)
    cluster_report = pd.read_csv(args.cluster_report)[["cluster_id","read_id"]]
    cluster_report.rename({'read_id': 'id'}, axis=1, inplace=True)
    
    # clustering discards single flnc reads;
    # subset only reads that are kept in for clustering to ease merging
    kept_flnc_reads = flnc_reads.loc[flnc_reads.id.isin(cluster_report["id"].values),]
    
    # merge cluster_report and flnc_reads 
    demuxed_cluster = pd.merge(cluster_report, kept_flnc_reads, how = "outer", on = "id")
    
    # write output
    demux_file_name = args.dir + "/" + args.output + ".clustered.demuxed.cluster_report.csv"
    print("Writing final files to", demux_file_name, "\n")
    demuxed_cluster.to_csv(demux_file_name,index=False)


def main():
    parser = argparse.ArgumentParser(description="Demultiplex merged cluster_report.csv using sample flnc read id")
    parser.add_argument("refine_dir", help='\t\tpath of refine directory containing <sample>_.flnc.report.csv')
    parser.add_argument('cluster_report', help='\t\tcsv file of cluster_report.csv after running IsoSeq3 cluster on all samples')
    parser.add_argument('-o','--output', required=False, help='\t\tPrefix for demuxed_cluster_report.csv; Default: prefix of cluster_report')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of cluster_report.csv.')
  
    args = parser.parse_args()
    if args.output is None:
        args.output = os.path.basename(args.cluster_report).split(".")[0] 
    
    if args.dir is None:
        args.dir =  os.path.dirname(args.cluster_report)
        
    demux_isoseq_cluster(args)
    print("All Done")
    

if __name__ == "__main__":
    main()
