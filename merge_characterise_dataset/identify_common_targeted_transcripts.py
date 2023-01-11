#!/usr/bin/env python3
# Szi Kay Leung (sl693@exeter.ac.uk)
# Parse through multiple files to find the number of transcripts commonly and uniquely detected by ONT and Iso-Seq 

import pandas as pd
from gtfparse import read_gtf
import csv
import sys
import argparse

## Require Target Gene List to convert ENSEMBL Id to Gene name 
TargetGenes_Ensembl_ID = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetGenesEnsembleId.txt"

def read_and_subset_target_genes(sqanti_input_file, abundance_file, dataset):
    TargetGene  = ["Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"]
    print("Reading:", sqanti_input_file)
    df = pd.read_csv(sqanti_input_file, sep = "\t",dtype={"bite": "string","polyA_motif":"string"}) 
    
    if dataset == "ONT":
        abundance = pd.read_csv(abundance_file, sep = "\t")
        Samples = ["S19","K24","L22","M21","O18","O23","O22","P19","T20","Q20","Q21","S18","S23","Q18","Q17","L18","Q23","T18"]
        abundance = abundance.loc[:, abundance.columns.isin(["annot_transcript_id"] + Samples)]
        abundance["FL"] = abundance.sum(axis=1)
        df.drop('FL', axis=1, inplace=True)
        df = df.merge(abundance, left_on='isoform', right_on='annot_transcript_id') 
    
    df = df[df['associated_gene'].isin(TargetGene)]
    return(df)


def read_input_files(args):
    IsoSeq = read_and_subset_target_genes(args.iso,"NA","Iso-Seq")
    ONT_Unfiltered = read_and_subset_target_genes(args.ont_1,args.a_ont,"ONT")
    ONT_Filtered =  read_and_subset_target_genes(args.ont_2,args.a_ont,"ONT")     
    Cuff_tmap = pd.read_csv(args.cuff, sep = "\t") 
    
    #ONT_Abundance = pd.read_csv(args.a_ont_1, sep = "\t")    
    #ONT_Filtered_Abundance = pd.read_csv(args.a_ont_2, sep = "\t")
    #ONT_UnFiltered_SQANTI_Abundance = pd.merge(ONT_Abundance,ONT_SQANTI, left_on="annot_transcript_id",right_on="isoform",how="left")
    #ONT_UnFiltered_SQANTI_Abundance.drop('isoform', axis=1, inplace=True) # to avoid downstream merge complications
    
    return IsoSeq, ONT_Unfiltered, ONT_Filtered, Cuff_tmap


def write_output_id(output_file, ID_List):
    textfile = open(output_file, "w")
    print("Writing output to:", output_file)
    for element in ID_List:
        textfile.write(element + "\n")
    textfile.close()


def identify_ID_and_stats(IsoSeq,Cuff_tmap,ONT_Unfiltered,ONT_Filtered,output_dir):
   
    ONT_Unfiltered_IDs = list(ONT_Unfiltered["isoform"])
    ONT_Filtered_IDs = list(ONT_Filtered["isoform"])
    
    Cuff_tmap_exact = Cuff_tmap[[True if i in ["="] else False for i in Cuff_tmap.class_code]]
    
    # list of PacBio transcripts that are identically matched to ONT transcriptome
    PB_detected = list(Cuff_tmap_exact.ref_id.unique())
    PB_detected_TargetGenes = list(set(PB_detected).intersection(IsoSeq["isoform"]))
    
    # filter for Target Genes 
    Cuff_tmap_exact = Cuff_tmap_exact[Cuff_tmap_exact['ref_id'].isin(PB_detected_TargetGenes)]
    
    # ONT Transcripts that are exact match to Iso-Seq (independent of Filtered or not)     
    ONT_ID_Match = list(Cuff_tmap_exact["qry_id"].unique())
    
    # ONT Filtered Transcripts that are exact match to Iso-Seq 
    ONT_Filtered_ID_Match = list(Cuff_tmap_exact.loc[Cuff_tmap_exact["qry_id"].isin(ONT_Filtered_IDs)].qry_id.unique())
    
    # Saving the reference of the standard output
    #sys.stdout = open(output_dir + "/" + "Compare_Stats.txt", 'w')
    print("***Total number***")
    print("Iso-Seq transcripts:", len(IsoSeq["isoform"].unique()))
    print("ONT unfiltered transcripts:", len(ONT_Unfiltered_IDs))
    print("ONT filtered transcripts:", len(ONT_Filtered_IDs))
    print("Matched Iso-Seq transcripts:", len(PB_detected_TargetGenes))
    print("Matched ONT transcripts:", len(ONT_ID_Match))
    print("Unique Iso-Seq transcripts:", len(set(IsoSeq["isoform"]) - set(PB_detected_TargetGenes)))
    print("Matched ONT filtered transcripts:", len(ONT_Filtered_ID_Match))
    #sys.stdout.close()
    
    # ONT Filtered Transcripts, Iso-Seq transcripts, 
    # Retaining all Iso-Seq transcripts (which would include the ones that are also in the unfiltered dataset)
    All_Retained = list(set(ONT_Filtered_IDs + list(IsoSeq["isoform"].unique())))
    ONT_All_Retained = list(set(ONT_Filtered_IDs + list(Cuff_tmap_exact["qry_id"])))
    print("All ONT Retained:", len(ONT_All_Retained))
    
    write_output_id(output_dir + "/ONT_MatchedIDs.txt", ONT_ID_Match)
    write_output_id(output_dir + "/ONT_Filtered_MatchedIDs.txt", ONT_Filtered_ID_Match)
    write_output_id(output_dir + "/IsoSeq_MatchedIDs.txt", PB_detected_TargetGenes)
    write_output_id(output_dir + "/IsoSeqONT_RetainedIDs.txt", ONT_All_Retained)
    write_output_id(output_dir + "/All_RetainedIDs.txt", All_Retained)
    Cuff_tmap_exact.to_csv(output_dir + "/CuffcompareExact.txt")
    Cuff_tmap_exact.to_csv(output_dir + "/ONT_All_Retained.txt")
    
    return All_Retained, ONT_All_Retained, Cuff_tmap_exact

    
    
def writeGTF(inGTF, file_path):
    
    '''
    #https://github.com/mpg-age-bioinformatics/AGEpy/blob/master/AGEpy/gtf.py
    Aim: write a GTF dataframe into a file
    :inGTF = GTF dataframe to be written. 
    :file_path = path/to/the/file.gtf
    '''
    
    cols=inGTF.columns.tolist()
    if len(cols) == 9:
        if 'attribute' in cols:
            df=inGTF
    else:
        df=inGTF[cols[:8]]
        df['attribute']=""
        for c in cols[8:]:
            if c == cols[len(cols)-1]:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'";'
            else:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'"; '
                
    df.to_csv(file_path, sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)


def subset_gtf(retainedID, sqanti_input_gtf, output_dir):

    All_Retained = pd.DataFrame({'Transcripts':retainedID})
    
    gffcompare_gtf = read_gtf(sqanti_input_gtf) 
    output_name = sqanti_input_gtf.split("/")[-1]
    output_name = output_dir + "/" + output_name.split("_")[0] + "_Final.gtf"
    
    NewId_Retained = list(gffcompare_gtf[gffcompare_gtf['transcript_id'].isin(All_Retained["Transcripts"])]["transcript_id"])
    FilteredGtf = gffcompare_gtf[gffcompare_gtf['transcript_id'].isin(NewId_Retained)]
    FilteredGtf = FilteredGtf[["seqname","source","feature","start","end","score","strand","frame","transcript_id","gene_id"]]
    
    print("Writing:", output_name)  
    writeGTF(FilteredGtf, output_name)
    
    return(FilteredGtf)
    

def merge_gtf(IsoSeq, Iso_SubsetGtf, ONT_SubsetGtf, Cuff_tmap_exact, output_dir):
    Cuff_tmap_exact["matched_id"] = Cuff_tmap_exact["ref_id"] + str("_") + Cuff_tmap_exact["qry_id"]
    Matched_dict = dict(zip(Cuff_tmap_exact["ref_id"],Cuff_tmap_exact["matched_id"]))
    IsoSeqFinal = pd.concat([Iso_SubsetGtf,ONT_SubsetGtf])
    IsoSeqFinal['transcript_id'].replace(Matched_dict, inplace=True)
    
    # remove ONT transcripts that are already matched with IsoSeq to avoid redundancies
    IsoSeqFinal = IsoSeqFinal[~IsoSeqFinal['transcript_id'].isin(Cuff_tmap_exact["qry_id"])] 
    
    
    print("Writing:", output_dir + "/IsoSeqONT_final.gtf")
    writeGTF(IsoSeqFinal, output_dir + "/IsoSeqONT_final.gtf")
    
    ## Convert gene id to gene name 
    IsoSeqGene_dict = dict(zip(["PB." + i.split(".",2)[1] for i in IsoSeq["isoform"]],IsoSeq["associated_gene"]))
    TargetGene = pd.read_csv(TargetGenes_Ensembl_ID, sep = "\t", header = None)
    TargetGene_dict = dict(zip(TargetGene[0],TargetGene[1]))
    IsoSeqFinal['gene_id'].replace(TargetGene_dict, inplace=True)
    IsoSeqFinal['gene_id'].replace(IsoSeqGene_dict, inplace=True)
    
    print("Writing:", output_dir + "/IsoSeqONT_final_genename.gtf")    
    writeGTF(IsoSeqFinal, output_dir + "/IsoSeqONT_final_genename.gtf")
    
    return(IsoSeqFinal)


def tabulate_abundance(IsoSeqFinal, ONT_Unfiltered, ONT_Filtered, IsoSeq, Cuff_tmap_exact,ONT_All_Retained, output_dir):
    
    IsoSeqFinal = IsoSeqFinal[IsoSeqFinal["feature"]=="transcript"].reset_index()
    #print(list(ONT_Unfiltered))
    
    ## ONT
    ONTcol = ["annot_transcript_id",'S19','K24','L22','M21','O18','O23','O22','P19','T20','Q20','Q21','S18','S23','Q18','Q17','L18','Q23','T18',"ONT_FL"]
    
    # Abundance of all other samples    
    # 1. ONT_unfiltered_retained = SQANTI classification with abundance of ONT transcripts retained in the final downstream
    # 2. ONT_unfiltered = SQANTI classification with abundance of all ONT unfiltered transcripts
    ONT_Unfiltered_retained = ONT_Unfiltered.loc[ONT_Unfiltered["isoform"].isin(ONT_All_Retained),] 
    ONT_Unfiltered_retained.to_csv(output_dir + "/ONT_retained_classification.csv")
    ONT_Unfiltered.to_csv(output_dir + "/ONT_Unfiltered_Abundance.csv", index=False)  
    
    #ONT_UnFiltered_SQANTI_Abundance["ONT_FL"] = ONT_UnFiltered_SQANTI_Abundance[Samples].sum(axis=1)
    ONT_Unfiltered.rename(columns={'FL': 'ONT_FL'}, inplace=True)
    ONT_Abundance = pd.merge(IsoSeqFinal[["transcript_id"]],ONT_Unfiltered[ONT_Unfiltered.columns.intersection(ONTcol)], left_on = "transcript_id", right_on = "annot_transcript_id", how = "left").dropna()
    ONT_Unfiltered.to_csv(output_dir + "/ONT_Unfiltered_Abundance.csv", index=False)   
    ONT_Unfiltered_retained[ONT_Unfiltered_retained.columns.intersection(ONTcol)].to_csv(output_dir + "/ONT_Retained_Abundance.csv", index=False)  
    #ONT_Abundance.to_csv(output_dir + "/ONT_Final_Abundance.csv") 
    
    
    ## IsoSeq
    IsoSeq_Samples = ['K17','K18','K19','K20','K21','K23','K24','L18','L22','M21','O18','O22','O23','P19','Q17','Q18','Q20','Q21','Q23','S18','S19','S23','T18','T20']
    IsoSeq["IsoSeq_FL"] = IsoSeq[["FL." + x for x in IsoSeq_Samples]].sum(axis=1)
    IsoSeq_Abundance = pd.merge(IsoSeqFinal[["transcript_id"]],IsoSeq[["isoform","IsoSeq_FL"]], left_on = "transcript_id", right_on = "isoform", how = "left").dropna() 
    IsoSeq_Abundance.to_csv(output_dir + "/IsoSeq_Abundance.csv")
    
    ## Merged         
    Merged_Abundance = pd.merge(Cuff_tmap_exact[["ref_id","qry_id"]],IsoSeq[["isoform","IsoSeq_FL"]], left_on="ref_id",right_on="isoform",how="left")
    Merged_Abundance = pd.merge(Merged_Abundance, ONT_Unfiltered[["annot_transcript_id","ONT_FL"]],left_on="qry_id",right_on="annot_transcript_id",how="left")
    Merged_Abundance["FL"] = Merged_Abundance[["IsoSeq_FL","ONT_FL"]].sum(axis=1)
    Merged_Abundance["isoform"] = Merged_Abundance["ref_id"] + str("_") + Merged_Abundance["qry_id"]
    Merged_Abundance.to_csv(output_dir + "/Matched_Abundance.csv")
    Final_Merged_Abundance = pd.merge(IsoSeqFinal[["transcript_id"]],Merged_Abundance[["isoform","FL"]], left_on = "transcript_id", right_on = "isoform", how = "left").dropna()
    
    # Concatenated
    ONT_Abundance = ONT_Abundance.rename({"ONT_FL":"FL"},axis=1)
    IsoSeq_Abundance = IsoSeq_Abundance.rename({"IsoSeq_FL":"FL"},axis=1)
    Concat_Abundance = pd.concat([ONT_Abundance,IsoSeq_Abundance,Final_Merged_Abundance])
    Concat_Abundance = Concat_Abundance[["transcript_id","FL"]]
    Concat_Abundance.rename({'transcript_id': 'id','FL': 'MergedFL'}, axis=1, inplace=True)
    Concat_Abundance["MergedFL"].astype(int)
    print("Writing:", output_dir + "/Final_Merged_Abundance.csv")
    Concat_Abundance.to_csv(output_dir + "/Final_Merged_Abundance.csv", index = False)   
    
    print("Total number of transcripts for downstream annotations:", len(Concat_Abundance))    
    return(Concat_Abundance)



def main():
    parser = argparse.ArgumentParser(description="Identifying Transcripts from ONT and Iso-Seq Targeted Transcriptome for annotation")
    parser.add_argument('--iso', "--isoseq_sqanti_class", help='\t\tIso-Seq SQANTI classification output file.')
    parser.add_argument('--iso_gtf', "--isoseq_sqanti_gtf", help='\t\tIso-Seq SQANTI classification gtf output file.')
    parser.add_argument('--ont_gtf',"--ont_unfiltered_sqanti_gtf", help='\t\tONT SQANTI classification output gtf file from TALON Unfiltered dataset.')
    parser.add_argument('--ont_1',"--ont_unfiltered_sqanti_class", help='\t\tONT SQANTI classification output file from TALON Unfiltered dataset.')
    parser.add_argument('--ont_2',"--ont_filtered_sqanti_class", help='\t\tONT SQANTI classification output file from TALON filtered dataset.')
    parser.add_argument('--a_ont', "--ont_unfiltered_abundance", help='\t\tONT TALON abundance file from Unfiltered dataset.')
    parser.add_argument('--cuff', "--cuff_reference_tmap", help='\t\tGffcompare cuff tmap output from IsoSeq as reference and ONT unfiltered dataset as annotation')
    parser.add_argument('--o_dir', "--output_dir", help='\t\tOutput path and name for list of ONT retained transcript IDs')

    args = parser.parse_args()
    IsoSeq, ONT_Unfiltered, ONT_Filtered, Cuff_tmap = read_input_files(args)
    All_Retained, ONT_All_Retained, Cuff_tmap_exact = identify_ID_and_stats(IsoSeq,Cuff_tmap,ONT_Unfiltered,ONT_Filtered,args.o_dir)
    Iso_SubsetGtf = subset_gtf(All_Retained, args.iso_gtf, args.o_dir)
    ONT_SubsetGtf = subset_gtf(All_Retained, args.ont_gtf, args.o_dir)
    IsoSeqFinal = merge_gtf(IsoSeq, Iso_SubsetGtf, ONT_SubsetGtf, Cuff_tmap_exact, args.o_dir)
    Final_abundance = tabulate_abundance(IsoSeqFinal, ONT_Unfiltered, ONT_Filtered, IsoSeq, Cuff_tmap_exact,ONT_All_Retained, args.o_dir)
    
    print("All Done")
    
    
if __name__ == "__main__":
    main()



