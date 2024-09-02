#!/usr/bin/env python

import os
import os.path
import subprocess 
import pandas as pd
import argparse
from argparse import Namespace
import distutils.spawn
import sys

old_stdout = sys.stdout

RSCRIPTPATH = distutils.spawn.find_executable('Rscript')
RSCRIPT_DIU = os.path.dirname(os.path.realpath(__file__)) + "/run_DIU_commandLine.R"
# include rownames in classification file

if os.system(RSCRIPTPATH + " --version")!=0:
    print("Rscript executable not found! Abort!", file=sys.stderr)
    sys.exit(-1)
    
    
def generate_gene_specific_files(args, gene):
  
    print("Generating gene-specific files:", gene)
    print("Writing output files to:", args.dir)
        
    # extract the associated isoform
    with open(args.classFile, 'r') as file:
        first_line = file.readline()  # Read the first line
        with open(args.dir + gene + '_classification.txt', 'w') as test1_file:  # Open in write mode to clear the file
            test1_file.write(first_line)

        with open(args.dir + gene + '_isoform.txt', 'w') as test2_file:  # Open in write mode to clear the file
            test2_file.write("")  # You can write an empty string or any initial content here

        for line in file.readlines()[1:]:  # Skip the first line
            associated_gene = line.split("\t")[first_line.split("\t").index("associated_gene")]
            isoform = line.split("\t")[0]
            if associated_gene == gene:
                with open(args.dir + gene + '_classification.txt', "a") as myfile:
                    myfile.write(line + "\n")
                with open(args.dir + gene + '_isoform.txt', "a") as myfile:
                    myfile.write(isoform + "\n")
                    
    # sanity check there are isoforms in the normalized expression file
    if os.stat(args.dir + gene + "_classification.txt").st_size == 0:
       print("ERROR: no isoforms subset in expression file")
       sys.exit()
  
    # note isoforms that are in the classfile might not be in the normalise expression file 
    print("Extract associated isoforms from normalised expression file")
    grepcmd = "grep -f {f} {n} > {o}".format(f=args.dir+gene+'_isoform.txt',n=args.normExp,o=args.dir+gene+'_normalised_expression.txt')
    sedcmd = "sed  -i '1i isoform,sample,normalised_counts' {f}".format(f=args.dir+gene+'_normalised_expression.txt')
    subprocess.run(grepcmd, shell=True)
    subprocess.run(sedcmd,shell=True)
    
    if os.stat(args.dir + gene + "_normalised_expression.txt").st_size == 0:
       print("ERROR: no isoforms subset in expression file")
       subprocess.run("rm {f}".format(f=args.dir+gene+'_normalised_expression.txt'), shell=True)
       sys.exit()
    
    

def run_DIU(args, gene):
    runRDIU = RSCRIPTPATH + " {s} -c {d}/{g}{c} -e {d}{g}{e} -t {t} -p {p} -f {f} &> {d}/{g}_r.log ".format(s=RSCRIPT_DIU,c='_classification.txt',e='_normalised_expression.txt',t=args.name,d=args.dir,g=gene,p=args.phenotype,f=args.factor)
    print(RSCRIPT_DIU)
    subprocess.run(runRDIU, shell=True)
    

def main():
    parser = argparse.ArgumentParser(description="Subset fasta and gtf from list of ID")
    parser.add_argument('classFile', help='\t\tSQANTI classification file (unfiltered)')
    parser.add_argument('normExp', help='\t\tNormalised expression matrix')
    parser.add_argument('-i', '--id_list', required=False, help='\t\tTxt file of list of ID to subset, no rownames and colnames')
    parser.add_argument('-I', '--ID_list', required=False, nargs='+', help='\t\tlist of ID to subset - direct on command line delimited by ,')
    parser.add_argument('-p', '--phenotype', help='\t\tInput phenotype file for downstream DIU analysis')
    parser.add_argument('-f', '--factor', help='\t\tInput factor file for downstream DIU analysis')
    parser.add_argument('-n', '--name', help='\t\tcommon character in expression columns to capture from classification file')
    parser.add_argument('-d','--dir', required=False, help='\t\tDirectory for output files. Default: Directory of input bed.')
    

    args = parser.parse_args()
    
    if args.dir is None:
       if not os.path.dirname(args.classFile):
         # If it's empty, set it to the current directory
          args.classFile = os.path.join(os.getcwd(), args.classFile)

       args.dir = os.path.dirname(args.classFile) + "/"
    else:
       args.dir = args.dir + "/"
  
    if args.id_list is not None:
        print("Reading in", args.id_list)
        id_list = pd.read_csv(args.id_list, sep = "\t", header = None)
        genes = list(id_list[0].values)
    elif args.ID_list is not None:
        genes = [i.split(",") for i in args.ID_list][0]
    else:
        print("Error: need to include --id_list or --ID_list")
        sys.exit()
        
    
    for g in genes:
        log_file = open(args.dir + g + "_py.log","w")
        sys.stdout = log_file
    
        generate_gene_specific_files(args, g)
        run_DIU(args,g)
        
        sys.stdout = old_stdout
        log_file.close()
    

if __name__ == "__main__":
    main()
