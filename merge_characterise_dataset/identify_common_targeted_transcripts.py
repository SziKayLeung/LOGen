#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Parse through gffcompare output to identify commonly matched isoforms between ONT and Iso-Seq dataset
Functions:
    read_input_files(args)
    group_isoforms(args, sq, Cuff_tmap)
    merge_gtf(args, allID, sq, Cuff_tmap_exact)
    classifier_dataset(id1, id2, id1name, id2name)
    tabulate_abundance(args, sq, allID, Cuff_tmap_exact)
Pre-requisites:
    Ensure the counts at the sample level are included in the sqanti classification files
    Run merge_talon_sqanti_forcounts.R for ONT dataset to include the counts
    Sample ID in the sqanti classification file and the IDs in the args.sample need to match
"""

# packages
import pandas as pd
from gtfparse import read_gtf
from datetime import date 
import csv
import sys
import argparse
import os


# LOGen modules
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sys.path.append(LOGEN + "miscellaneous")
sys.path.append(LOGEN + "merge_characterise_dataset")
sys.path.append(LOGEN + "compare_datasets")
import read_write as rw
import subset_fasta_gtf as sub
import subset_targetgenes_classfiles as sub_target
import dataset_identifer as di

# turn off warnings in this script (relating to loc)
pd.options.mode.chained_assignment = None  # default='warn'


"""
read input files <sqanti classification files> and save df as dictionary(pbid:sequence)
read gffcompare tmap output file and filter by target genes using ensembl ID
:param args: args.tgenes_ens, args.iso, args.ont_unfil, args.ont_fil, args.cuff
:returns sq: dictionary of sqanti classification files 
:returns Cuff_tmap: gffcompare tmap output file
"""
def read_input_files(args):
    
    # read sqanti classification files and subset by associated target genes
    # save classification files into dictionary
    sq = {'iso' : sub_target.read_and_subset_target_genes(args.tgenes_ens, args.iso),
        'ont_unfilter' : sub_target.read_and_subset_target_genes(args.tgenes_ens, args.ont_unfil),
        'ont_filter' : sub_target.read_and_subset_target_genes(args.tgenes_ens, args.ont_fil)}
    
    # extract only target gene (using ensembl ID) in gffcompare output
    TargetGene_file = pd.read_csv(args.tgenes_ens, sep = "\t", header = None)
    TargetGene_Ensembl = list(TargetGene_file[[0]].values.ravel())
    
    # read in gffcommpare tmap output
    # however do not filter by target genes using ref_gene_id, due to misannotations from TALON
    # resulting in over-filtering of isoforms associated to target gene
    Cuff_tmap = pd.read_csv(args.cuff, sep = "\t") 
            
    return sq, Cuff_tmap
    

"""
identify and group isoforms between Iso-Seq and ONT dataset using gffcompare output
matched isoforms defined as "=" in gffcompare output (tmap file)
:param args: args.o_dir
:returns allID: dictionary of groups of isoforms
:returns Cuff_tmap_exact: gffcompare tmap output file with exact match ("=") of isoforms
"""
def group_isoforms(args, sq, Cuff_tmap):
    
    # only extracting the isoforms with exact match "=" from gffcompare output
    Cuff_tmap_exact = Cuff_tmap[[True if i in ["="] else False for i in Cuff_tmap.class_code]]
    
    # create a column with new matched ID for matched transcripts
    Cuff_tmap_exact.loc[:,"matched_id"] = Cuff_tmap_exact["ref_id"] + "_" + Cuff_tmap_exact["qry_id"]
    
    # group and generate list of isoforms
    '''
    i) ID = list of isoforms from input file 
    Iso-Seq & ONT sqanti classification file
    Gffcompare tmap output file 
    ONT_unfilter_unique: ONT isoforms discarded using TALON filtering
    iso_match and ont_match: matching Iso-Seq and ONT isoforms associated with target genes
    hence use "intersection" (i.e common isoforms) with sqanti classification files that have already been subsetted with target genes
    NB: Iso-Seq gtf used as reference when running gffcompare, 
        therefore ref_id = Iso-Seq isoforms, qry_id = ONT isoforms
        ONT unfiltered gtf used for gffcompare
    
    ii) matchID = list of isoforms that are matched based on gffcompare output
    iso_match = Iso-Seq isoforms found in ONT dataset (unfiltered)
    ont_filter_match = ONT isoforms in ONT filtered dataset also found in Iso-Seq dataset
    ont_unfilter_match = ONT isoforms unique in ONT unfiltered dataset also found in Iso-Seq dataset 
    (i.e. otherwise would have discarded)
        
    iii) unmatchID = list of isoforms that not matched based on gffcompare output
    iso_unmatch = Iso-Seq isoforms not found in ONT dataset (unfiltered)
    ont_unfilter_unmatch = ONT isoforms discarded under TALON filtering and not found in Iso-Seq dataset
    ont_filter_unmatch = ONT isoforms in ONT filtered dataset not found in Iso-Seq dataset
    
    iv) finalID = list of isoforms from merging of two datasets
    final_ont = all ONT transcripts kept: 
                ONT filtered dataset & unique ONT isoforms not retained from filtering but matched in Iso-Seq dataset
    merge_iso = merged list of isoforms (Iso-Seq+ONT) using Iso-Seq ID from matched isoforms (to avoid redundancy)
                all Iso-Seq isoforms + ONT filtered but not matched isoforms
    merge_ont = merged list of isoforms (Iso-Seq+ONT) using ONT ID from matched isoforms (to avoid redundancy)
                unmatched Iso-Seq isoforms (unique Iso-Seq isoforms) + matched ONT isoforms + ONT filtered but not matched isoforms
    NB: len(merge_iso) == len(merge_ont)
    '''
    ID = {'iso' : list(sq["iso"]["isoform"]),
          'ont_unfilter' : list(sq["ont_unfilter"]["isoform"]),
          'ont_filter': list(sq["ont_filter"]["isoform"]),
          'ont_unfilter_unique' : list(set(sq["ont_unfilter"]["isoform"]) - set(sq["ont_filter"]["isoform"])),
          'iso_match' : list(set(Cuff_tmap_exact["ref_id"]).intersection(sq["iso"]["isoform"])),
          'ont_match': list(set(Cuff_tmap_exact["qry_id"]).intersection(sq["ont_unfilter"]["isoform"]))}

    matchID = {'iso_match' : ID["iso_match"],
        'ont_filter_match': list(set(ID["ont_match"]).intersection(ID["ont_filter"])),
        'ont_unfilter_match': list(set(ID["ont_match"]) - set(ID["ont_match"]).intersection(ID["ont_filter"]))}

    unmatchID = {'iso_unmatch' : list(set(ID["iso"]) - set(matchID["iso_match"])),
        'ont_unfilter_unmatch' : list(set(ID["ont_unfilter_unique"]) - set(matchID["ont_unfilter_match"])),
        'ont_filter_unmatch' : list(set(ID["ont_filter"]) - set(matchID["ont_filter_match"]))}

    finalID = {'final_ont' : list(matchID["ont_unfilter_match"] + ID["ont_filter"]),
        'merge_iso' : list(ID["iso"] + unmatchID["ont_filter_unmatch"]),
        'merge_ont' : list(unmatchID["iso_unmatch"] + ID["ont_match"] + unmatchID["ont_filter_unmatch"])}

    # combine all ID into one dictionary
    allID = {**ID, **matchID, **unmatchID, **finalID}

    # loop through combined dictionary <name = category; iso = list of isoforms>
    # extract the number of isoforms under each category; saved under nums
    # write the isoform ID into txt files
    nums = {}
    for name, iso in allID.items():
        nums[name] = str(len(iso))
        rw.write_lst(args.o_dir + "/isoform_id/" + name + ".txt", iso)
    nums
    
    # write nums to output file
    num_df = pd.DataFrame(nums.items(), columns=['key', 'num_isoforms'])
    num_df.to_csv(args.o_dir + "/num_isoforms.txt", index=False)
    
    return allID, Cuff_tmap_exact

    
"""
generate a finalised merged gtf based, removing redundant ONT matched isoforms
:param args: args.o_dir
:param allID: dictionary of groups of isoforms
:param sq: dictionary of read classification files
:param Cuff_tmap_exact: gffcompare tmap output file with exact match ("=") of isoforms
:returns Nothing
"""
def merge_gtf(args, allID, sq, Cuff_tmap_exact):
    
    '''
    1. generate Iso-Seq retained gtf and ONT retained gtf
    2. merge to create a unified retained gtf
    3. replace Iso-Seq ID of matched isoforms in unified gtf with matched ID containing Iso-seq and ONT ID
    4. remove redundant ONT matched isoforms in unified gtf
    5. replace gene_id column with gene names 
    
    finalised gtf consisting of isoforms from
    Iso-Seq/ONT matched: <Iso-Seq-ID>_<ONT_ID>
    Iso-Seq unmatched: <Iso-Seq-ID>
    ONT filtered unmatched: <ONT_ID>
        
    NB: gtf takes Iso-Seq gtf as reference for matched isoforms
    matched ID comes from Cuff_tmap_exact
    '''
    
    # generate dataset specific gtf 
    # Iso_retained_gtf = all Iso-Seq isoforms
    # ONT_retained_gtf = all ONT filtered isoforms + ONT discarded but matched isoforms
    Iso_retained_gtf = sub.subset_gtf(allID["iso"], args.iso_gtf, "retained", dir=args.o_dir)
    ONT_retained_gtf = sub.subset_gtf(allID["final_ont"], args.ont_gtf, "retained", dir=args.o_dir)
    
    # merge Iso-Seq and ONT retained gtf
    # NB: redundant isoforms at this point, whereby matching isoforms appear twice with Iso-Seq and ONT IDs
    All_retained_gtf = pd.concat([Iso_retained_gtf, ONT_retained_gtf])
    
    # create a dictionary: key = Iso-Seq Matched ID (ref_id), value = matched ID
    # no need to create a dictionary for ONT Matched ID
    # as ONT matched transcripts are removed to avoid redundancies
    # thereby using the PacBio transcript as reference for matched IDs
    Matched_dict = dict(zip(Cuff_tmap_exact["ref_id"],Cuff_tmap_exact["matched_id"]))
    
    # replace original Iso-Seq transcript ID with matched ID 
    All_retained_gtf['transcript_id'].replace(Matched_dict, inplace=True)
    
    # remove ONT transcripts that are already matched with IsoSeq to avoid redundancies
    # qry_id = ONT isoforms that are matching given ONT gtf is used as query gtf when running gffcompare
    All_retained_gtf = All_retained_gtf[~All_retained_gtf['transcript_id'].isin(Cuff_tmap_exact["qry_id"])] 
    
    print("Writing:", args.o_dir + "/IsoSeqONT_final.gtf")
    rw.writeGTF(All_retained_gtf, args.o_dir + "/IsoSeqONT_final.gtf")
    
    # QC 
    if not len(set(list(filter(lambda x:'_' in x, All_retained_gtf['transcript_id'])))) == len(allID["iso_match"]):
        print("Matching isoforms missing in gtf")
    if len(set(allID["iso_unmatch"]) - set(All_retained_gtf['transcript_id'])) > 0:
        print("Missing unmatched Iso-Seq isoforms in gtf")
    if len(set(allID["ont_filter_unmatch"]) - set(All_retained_gtf['transcript_id'])) > 0:
        print("Missing filtered by unmatched ONT isoforms in gtf")
    
    # Convert gene id to gene name 
    # gtf gene_name currently either <PB.X> if PacBio dataset of <ENS..> if ONT dataset
    # create IsoSeqGene_dict and replace gene_ID in gtf
    IsoSeqGene_dict = dict(zip(["PB." + i.split(".",2)[1] for i in sq["iso"]["isoform"]],sq["iso"]["associated_gene"]))
    All_retained_gtf['gene_id'].replace(IsoSeqGene_dict, inplace=True)
    
    # create TargetGene_dict with Ensembl ID and replace gene_ID in gtf
    TargetGene = pd.read_csv(args.tgenes_ens, sep = "\t", header = None)
    TargetGene_dict = dict(zip(TargetGene[0],TargetGene[1]))
    All_retained_gtf['gene_id'].replace(TargetGene_dict, inplace=True)
    
    print("Writing:", args.o_dir + "/IsoSeqONT_final_genename.gtf")    
    rw.writeGTF(All_retained_gtf, args.o_dir + "/IsoSeqONT_final_genename.gtf")


"""
tabulate abundance across Iso-Seq and ONT dataset (noting matched isoforms)
matched isoforms defined as "=" in gffcompare output (tmap file)
:param args: args.sample, args.o_dir
:param sq: dictionary of sqanti classification files read (Iso-Seq, ONT unfiltered)
:param allID : dictionary of grouped isoforms 
:param  Cuff_tmap_exact: df of matched isoforms from gffcompare output
:returns Nothing
"""
def tabulate_abundance(args, sq, allID, Cuff_tmap_exact):
    
    '''
    1. extract the sample ID from each dataset using input args.sample file
    2. create a final ONT classification file of the ONT isoforms retained (filtered + matched unfiltered)
    3. calculate sum and median of FL reads across samples per isoform in Iso-Seq and final ONT dataset; 
       store as ONT_sum_FL, ONT_med_FL, IsoSeq_sum_FL, IsoSeq_med_FL
    4. create a df with Iso-Seq isoforms and Cuff_tmap_exact to identify which Iso-Seq isoforms are matched
    5. add ONT isoforms to df to identify which ONT isoforms are matched
    6. create a "dataset" column in df using classifier_dataset()
    7. include abundance at a sample level, renaming sample columns with "Iso-Seq" and "ONT" in front to differentiate
    '''
    
    # read in file with sample IDs
    samples = pd.read_csv(args.sample)
    
    # remove NA from list and store in dictionary
    samples = {'ont' : [x for x in samples['ONT'].values if str(x) != 'nan'],
               'iso' : [x for x in samples['IsoSeq'].values if str(x) != 'nan']}
    
    # create final ONT sqanti classification file by subsetting the ONT isoforms that are retained 
    sq['final_ont'] = sq['ont_unfilter'].loc[sq['ont_unfilter']['isoform'].isin(allID['final_ont']),] 

    # create two new columns of the sum and median FL reads across samples
    sq['final_ont'].loc[:,"ONT_sum_FL"] =  sq['final_ont'][samples['ont']].sum(axis=1)
    sq['final_ont'].loc[:,"ONT_med_FL"] =  sq['final_ont'][samples['ont']].median(axis=1)
    
    # NB: Iso-Seq sqanti classification file might have "FL.X" columns, where X is a sample ID
    # therefore if "FL." is noted in the column names, then paste "FL." in the samples 
    # before finding the sum and median
    if len(sq["iso"].filter(like='FL.').columns) > 0:
        sq['iso'].loc[:,"IsoSeq_sum_FL"] =  sq['iso'][["FL." + x for x in samples['iso']]].sum(axis=1)
        sq['iso'].loc[:,"IsoSeq_med_FL"] =  sq['iso'][["FL." + x for x in samples['iso']]].median(axis=1)
        # replace isoseq dataset of samples columns with "FL.X" to "X"
        sq['iso'].columns = sq['iso'].columns.str.replace('FL.',"")
    else:
        sq['iso'].loc[:,"IsoSeq_sum_FL"] =  sq['iso'][samples['iso']].sum(axis=1)
        sq['iso'].loc[:,"IsoSeq_med_FL"] =  sq['iso'][samples['iso']].median(axis=1)
    
    # merge abundance first for Iso-Seq with Cuff_tmap_exact 
    # need to merge to determine which Iso-Seq isoforms are matched and found in the ONT dataset
    # merge on "left" to ensure all the Iso-Seq isoforms detected are also captured
    # right_on = "ref_id" as Iso-Seq gtf used to run gffcompare (therefore ref_id refers to Iso-Seq isoforms)
    iso_abundance = pd.merge(sq['iso'][["isoform","IsoSeq_sum_FL","IsoSeq_med_FL"]], Cuff_tmap_exact[["ref_id","qry_id","matched_id"]], 
                         left_on = "isoform", right_on="ref_id", how = "left")

    # merge abundance second to include ONT isoforms in previous abundance file
    # merge isoforms from ONT final classification file with the qry_id column 
    # qry_id column kept from previously merging with Cuff_tmap_exact
    # merge on "outer" to keep all the ONT isoforms that are not matched (i.e does not have a qry_id in cuff_tmap_exact)
    # and to keep all Iso-Seq isoforms in the iso_abundance
    abundance = pd.merge(sq['final_ont'][["isoform","ONT_sum_FL","ONT_med_FL"]], iso_abundance,
                     left_on = "isoform", right_on="qry_id", how = "outer")
    
    # rename columns and drop unnecessary columns
    abundance.columns  = ["ont_isoform","ONT_sum_FL","ONT_med_FL","isoseq_isoform","IsoSeq_sum_FL","IsoSeq_med_FL","ref_id","qry_id","matched_id"]
    abundance.drop(['ref_id', 'qry_id'], axis=1)
    
    # running through each row, determine if isoform is "Both", "ONT", "Iso-Seq" 
    abundance["dataset"] = abundance.apply(
        lambda row: di.classifier_dataset(row['ont_isoform'], row['isoseq_isoform'],"ONT","Iso-Seq"),
         axis=1)
    
    # running through each row, create a unionised column depending on the dataset
    abundance["union_isoform"] = abundance.apply(
        lambda row: di.unionise_id(row["dataset"], row["isoseq_isoform"], row["ont_isoform"], row['matched_id']),
        axis=1)
     
    # re-arrange columns
    abundance = abundance[["dataset","union_isoform","isoseq_isoform","ont_isoform","IsoSeq_sum_FL","ONT_sum_FL","IsoSeq_med_FL","ONT_med_FL","matched_id"]]
    
    # include abundance at a sample level
    # generate a dictionary of the sqanti classification df with only the abundance columns + "isoform" column
    s_abundance = {
    'iso' : sq["iso"][["isoform"] + samples["iso"]],
    'ont' : sq["final_ont"][["isoform"] + samples["ont"]]
    }
    
    # add "Iso-Seq" and "ONT" in front of the abundance column names to differentiate where the samples are from when merging later
    s_abundance["iso"].columns = ["isoseq_" + x for x in s_abundance["iso"].columns]
    s_abundance["ont"].columns = ["ont_" + x for x in s_abundance["ont"].columns]
    
    # merge the abundance at the sample level with the abundance df originally generated 
    # 1st merge the original abundance file with the Iso-Seq abundance at sample level; 
    # merge on "left" to keep all the isoforms in the original abundance file
    # 2nd merge the now updated Iso-Seq abundance with the ONT abundance at sample level;
    # merge on "left" to keep all the isoforms that were in the original abundance file, and the Iso-Seq isoforms
    s_abundance["iso"] = pd.merge(abundance, s_abundance["iso"], on = "isoseq_isoform", how = "left")
    final_abundance = pd.merge(s_abundance["iso"], s_abundance["ont"], on = "ont_isoform", how = "left")
    
    # QC: check the final abundance files has all the isoforms kept from Iso-Seq dataset and final ONT dataset
    # remove "na" from final_abundance, as there would be NAs from isoforms that are only unique to each dataset
    if len(set(list(final_abundance["ont_isoform"].dropna())) - set(list(sq["final_ont"]["isoform"]))) > 0:
        print("ONT isoforms not retained in abundance file")
        sys.exit()
        
    if len(set(list(final_abundance["isoseq_isoform"].dropna())) - set(list(sq["iso"]["isoform"]))) > 0:
        print("Iso-Seq isoforms not retained in abundance file")
        sys.exit()
        
    print("Writing:", args.o_dir + "/Final_Merged_Abundance.csv")
    final_abundance.to_csv(args.o_dir + "/Final_Merged_Abundance.csv", index = False)   


def main():
    parser = argparse.ArgumentParser(description="Identifying Transcripts from ONT and Iso-Seq Targeted Transcriptome for annotation")
    parser.add_argument('--iso', "--isoseq_sqanti_class", help='\t\tIso-Seq SQANTI classification output file.')
    parser.add_argument('--iso_gtf', "--isoseq_sqanti_gtf", help='\t\tIso-Seq SQANTI classification gtf output file.')
    parser.add_argument('--ont_gtf',"--ont_unfiltered_sqanti_gtf", help='\t\tONT SQANTI classification output gtf file from TALON Unfiltered dataset.')
    parser.add_argument('--ont_unfil',"--ont_unfiltered_sqanti_class", help='\t\tONT SQANTI classification output file from TALON Unfiltered dataset.')
    parser.add_argument('--ont_fil',"--ont_filtered_sqanti_class", help='\t\tONT SQANTI classification output file from TALON filtered dataset.')
    parser.add_argument('--cuff', "--cuff_reference_tmap", help='\t\tGffcompare cuff tmap output from IsoSeq as reference and ONT unfiltered dataset as annotation')
    parser.add_argument('--o_dir', "--output_dir", help='\t\tOutput path and name for list of ONT retained transcript IDs')
    parser.add_argument('--tgenes_ens', help='\t\tTxt file containing target genes; column 1 - ENSEMBL ID, column 2 - gene name')
    parser.add_argument('--sample', help='\t\tTxt file containing sampe IDs; column 1 - Iso-Seq, column 2 - ONT')

    args = parser.parse_args()
    
    # generate params file
    args.doc = os.path.join(os.path.abspath(args.o_dir) +"/params.txt")
    tgenes = ', '.join(pd.read_csv(args.tgenes_ens,sep="\t",header=None)[1].values)
    print("Write arguments to {0}...".format(args.doc, file=sys.stdout))
    with open(args.doc, 'w') as f:
        f.write("Date processed\t" + str(date.today()) + "\n")
        f.write("Input sqanti classification files:\n")
        f.write("Iso-Seq\t" + args.iso + "\n")
        f.write("ONT unfiltered\t" + args.ont_unfil + "\n")
        f.write("ONT filtered\t" + args.ont_fil + "\n")
        f.write("Input gffcompare\t" + args.cuff + "\n")
        f.write("Target genes\t" + tgenes + "\n")
        
    # make id directory 
    if not os.path.exists(args.o_dir + "/isoform_id/"):
        os.mkdir(args.o_dir + "/isoform_id/")
    
    # read input files
    sq, Cuff_tmap = read_input_files(args)
    
    # generate a dictionary of IDs based from grouping of isoforms (if matched etc)
    allID, Cuff_tmap_exact = group_isoforms(args, sq, Cuff_tmap)
    
    # generate gtf and abundance file
    merge_gtf(args, allID, sq, Cuff_tmap_exact)
    tabulate_abundance(args, sq, allID, Cuff_tmap_exact)
    
    print("All Done")
    
    
if __name__ == "__main__":
    main()
