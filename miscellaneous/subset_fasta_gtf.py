#!/usr/bin/env python

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Subset gtf based on transcript ID and read output
Functions:
    subset_gtf
"""

# packages
import os 
from gtfparse import read_gtf
import read_write as rw

"""
Subset gtf based on transcript ID 
:param retainedID: list of transcripts  
:param input_gtf: path of input gtf for subsetting
:param output name: string of name to be included in the output
:param *args[0]: optional argument, if present then generate the output file to this directory
:returns output_dict: dictionary of read name as key and sequence as value
"""
def subset_gtf(retainedID, input_gtf, name, **kwargs):
    
    print("Reading gtf:", input_gtf)
    gtf = read_gtf(input_gtf) 
    
    # variables
    basename = os.path.basename(input_gtf)
    outname = basename.replace(".gtf", "_" + name +".gtf")
    optdir = kwargs.get('dir', None)
    if optdir is not None:
        outdir = optdir
    else:
        outdir = os.path.dirname(input_gtf)
 
    subset_gtf = gtf[gtf['transcript_id'].isin(retainedID)]
    subset_gtf = subset_gtf[["seqname","source","feature","start","end","score","strand","frame","transcript_id","gene_id"]]
    
    print("Number of unique transcripts subsetted:", len(set(subset_gtf["transcript_id"].values)))
    print("Writing:", outdir + "/" + outname)  
    rw.writeGTF(subset_gtf, outdir + "/" + outname)
    
    return(subset_gtf)
