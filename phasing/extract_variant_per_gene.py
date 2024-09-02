import pandas as pd
import subprocess 
import os
import argparse
import sys
import shutil


def run_plink_per_gene(args):
    neuroChip = pd.read_csv(args.neuroChipFile, sep = "\t")
    snp = neuroChip.loc[neuroChip["ANNO_Gene.refGene"] == args.gene,"NeuroChip_variant_location_hg19"]
    snp.to_csv(args.odir + "snp.txt", sep='\t', encoding='utf-8', index=False, header=False)

    print("run plink")
    plinkExtract = "plink --bfile {i}/{b} --extract {n} --make-bed --out {d}{o}".format(b=args.b,n=args.odir+"snp.txt",o=args.gene,d=args.odir,i=args.dir)
    plinkRecode = "plink --bfile {d}/{b} --recode A --out {d}{o}".format(b=args.gene,o=args.gene+"_Variant",d=args.odir)
    plinkRecodeBasic = "plink --bfile {d}/{b} --recode --out {d}{o}".format(b=args.gene,o=args.gene+"_Variant",d=args.odir)
    subprocess.run(plinkExtract, shell=True)
    
    if os.path.isfile(args.odir + args.gene + ".bed"):
        print("Variant recorded")
        subprocess.run(plinkRecode, shell=True)
        subprocess.run(plinkRecodeBasic, shell=True)


def match_phenotype(args):
    if os.path.isfile(args.odir + args.gene + "_Variant.raw"):
        raw = pd.read_csv(args.odir + args.gene + "_Variant.raw", sep=' ')
        raw = raw.drop(columns=['FID', 'PAT','MAT','SEX','PHENOTYPE'])
        phenotype = pd.read_csv(args.phenotype)
        phenotype = phenotype[["BBNId","Brain_ID","DNA_ID"]]
        raw = phenotype.merge(raw, left_on = "DNA_ID", right_on = "IID")
        raw.to_csv(args.odir + args.gene + "_Variant_Phenotype.txt", sep='\t', encoding='utf-8', index=False, header=True)

    

def main():
    parser = argparse.ArgumentParser(description="Subset fasta and gtf from list of ID")
    parser.add_argument('neuroChipFile', help='\t\tInput neuroChip file containing variants')
    parser.add_argument('-g', '--gene', required=False, help='\t\tGene to extract snps')
    parser.add_argument('-b', help='\t\tplink prefix files (--bfile)')
    parser.add_argument('-p', '--phenotype', required=False, help='\t\tPhenotype data')
    parser.add_argument('-d','--dir', help='\t\tDirectory containing plink files, also generates output files.')

    args = parser.parse_args()
    args.odir = args.dir + "/" + args.gene + "/"
    
    if not os.path.exists(args.odir):
        os.makedirs(args.odir)
    else:
        print("directory already exists, replacing files.".format(args.odir), file=sys.stderr)
        shutil.rmtree(args.odir) 
        os.makedirs(args.odir)
    
    run_plink_per_gene(args)
    
    if args.phenotype is not None:
        match_phenotype(args)
    
if __name__ == "__main__":
    main()
