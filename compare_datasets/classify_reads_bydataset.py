#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Date: March 2023
Aim: Classify the isoform in the abundance file (each row) using the column names (samples), and create a dataset column

--- Input: 
1. Abundance ==> demux_FL_count.csv: row = isoform ID, column = sample ID
2. Phenotype ==> phenotye.csv: <sample> <group> 
# sample in phenotype.csv needs to match column names in abundance file

--- Output:
Phenotype file with "dataset" column

"""

# packages
import argparse 
import pandas as pd
import numpy as np
import os


"""
create_dataset_column
create a dataset column with the dataset names if read is detected in said dataset
:param args: args.abundance, args.phenotype
:returns final = df of args.abundance with additional columns: <datasetX>_sum_FL, and dataset 
"""
def create_dataset_column(args):
    '''
    1. read in abundance and phenotype file
        ==> abundance:
        id,sample1,sample2,sample3,sample4,..
        PB.1.1,0,2,3,0,..
        PB.2.1,0,0,1,0,..

        ==> phenotype:
        sample,group
        sample1,control
        sample2,control
        sample3,AD,
        sample4,intermediate
    
    2. datawrangle abundance file from wide to long, and merge with phenotype file
        ==> abLong:
        sample,isoform,FL,group
        sample1,PB.1.1,0,control
        sample2,PB.1.1,2,control
        sample3,PB.1.1,3,AD
        sample4,PB.1.1,0,intermediate
        sample1,PB.2.1,0,control
        sample2,PB.2.1,0,control
        sample3,PB.2.1,1,AD
        sample4,PB.2.1,0,intermediate
        
    3. sum FL reads across isoform and group
        ==> abTally:
        groupisoform,AD,control,intermediate
        PB.1.1,3,2,0
        PB.2.1,1,0,0
        
    4. create a dataset column based on the presence of reads in each column 
        ==> abTally:
        groupisoform,AD,control,intermediate,dataset
        PB.1.1,3,2,0,[AD,control]
        PB.2.1,1,0,0,[AD]
        e.g. PB.1.1 only detected in AD and control, whereas PB.2.1 detected only in AD
        Detection based if reads != 0
    
    5. Replace dataset name as "All" if reads are detected across all the groups, by
       matching the dataset name with the column name merged (columns = dataset names)
    
    6. Add "sum_FL" to column names if column names != "dataset"
    7. merge abTally dataset to original abundance file using index 
        note: PB. stripped from isoform column for faster merging       
    '''
    # 1. read in arguments
    print("Reading in files...")
    abundance = pd.read_csv(args.abundance)
    phenotype = pd.read_csv(args.phenotype)
    
    # 2. abundance df: wide to long
    abLong = pd.melt(abundance, id_vars='id')
    abLong.columns = ["isoform","sample","FL"]
    abLong = abLong.set_index('sample').join(phenotype.set_index('sample'), how='left') 
    abLong = abLong.reset_index()
    
    # 3. sum FL reads across phenotype group and isoform
    print("Sum reads...")
    abTally = abLong.groupby(['group','isoform']).sum("FL")
    abTally = abTally.reset_index()
    abTally = abTally.pivot(index='isoform', columns='group', values='FL')
    
    # 4. create a dataset column based on the presence of reads
    # alldataset = string of all dataset names concatenated and delimited by ","; for downstream matching
    alldataset = ','.join(abTally.columns)
    abTally["dataset"]=[','.join(filter(None, i)) for i in np.where(abTally.values==0,'',abTally.columns)]
    
    # 5. Search dataset column for alldataset value, and replace with "All"
    abTally["dataset"]=np.where(abTally['dataset']==alldataset, "All", abTally['dataset'])
    
    # 6. Add "sum_FL" to column names    
    abTally.columns = [i + "_sum_FL" if  i != "dataset" else i for i in abTally.columns]
    
    # 7. merge by index 
    print("Final merging...")
    abTally.index = abTally.index.str.replace("PB.","").astype(float)
    abundance.index = abundance["id"].str.replace("PB.","").astype(float)
    final = pd.merge(abundance, abTally, left_index=True, right_index=True)
    
    return(final)
    

def create_dataset_column(abLong):
    # 3. sum FL reads across phenotype group and isoform
    abTally = abLong.groupby(['group','isoform']).sum("FL")
    abTally = abTally.reset_index()
    abTally = abTally.pivot(index='isoform', columns='group', values='FL')
    
    # 4. create a dataset column based on the presence of reads
    # alldataset = string of all dataset names concatenated and delimited by ","; for downstream matching
    alldataset = ','.join(abTally.columns)
    abTally["dataset"]=[','.join(filter(None, i)) for i in np.where(abTally.values==0,'',abTally.columns)]
    
    # 5. Search dataset column for alldataset value, and replace with "All"
    abTally["dataset"]=np.where(abTally['dataset']==alldataset, "All", abTally['dataset'])
    
    # 6. Add "sum_FL" to column names    
    abTally.columns = [i + "_sum_FL" if  i != "dataset" else i for i in abTally.columns]
    
    return(abTally)
    
    
def read_line_byline(args):
  
    phenotype = pd.read_csv(args.phenotype)
    output = open(args.dir + "/" + args.output + "_classified.csv","w")
    
    with open(args.abundance) as file:
      cols=file.readline()
      cols=cols.replace("\n","")
      for index, line in enumerate(file):
          line = line.replace("\n","")
          df = pd.DataFrame({"sample" : cols.split(","), "FL" : line.split(",")})
          iso = df.iloc[0,1]
          df.drop(df.index[0],inplace=True)
          df = df.set_index('sample').join(phenotype.set_index('sample'), how='left') 
          df = df.reset_index()
          df["isoform"] = iso
          df["FL"]=pd.to_numeric(df["FL"], downcast='float')
          df2 = create_dataset_column(df)
          if index == 1:
              newcols = ",".join(df2.columns).replace("\n","")
              output.write(cols + "," + newcols + "\n")
          else:
              newline = ",".join([str(i) for i in df2.values[:1][0]])
              output.write(line + "," + newline + "\n")
      
    output.close()


def main():
    parser = argparse.ArgumentParser(description="Colour Transcripts from ONT and Iso-Seq Targeted Transcriptome by abundance and coding potential")
    parser.add_argument('-a',"--abundance", help='\t\tAbundance file of FL reads from final merged dataset')
    parser.add_argument('-p',"--phenotype", help='\t\tPhenotype file')
    parser.add_argument('-o',"--output", required=False,help='\t\tOutput coloured bed12 file from the final merged dataset')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of input bed.')
    
    args = parser.parse_args() 
    
    #---- QC arguments
    
    if args.dir is None:
        args.dir = os.path.dirname(args.abundance) 
    
    if args.output is None:
        args.output = os.path.basename(args.abundance).split(".")[0]

    print("Writing output files to:", args.dir + "/" + args.output + "_classified.csv")

    #---- run  
    abundance = create_dataset_column(args)
    abundance.to_csv(args.dir + "/" + args.output + "_classified.csv",index=False)
    #read_line_byline(args)
        
if __name__ == "__main__":
    main()
