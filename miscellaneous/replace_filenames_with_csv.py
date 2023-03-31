#!/usr/bin/env python

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Aim: Replace (or copy and replace) multiple file names in a directory using reference csv file
Input:
    file = csv file with 2 columns with header names old_name,new_name
    make sure old_name column is unique 
Notes:
    function loop through each file in input directory to cross-reference
    thus time consuming and looping if many files/hidden files etc
    ignores extension and replaces/copies all filenames that match
    
Functions:
    rename_file(args)
"""

# packages
import os
import pandas as pd
import argparse
import sys
import shutil

def rename_func(args, indexdict, fullname, filename):
        
    if args.ext:
        basename = filename.split("." + args.ext,)[0]
        extension = "." + args.ext
    else:
        basename, extension = os.path.splitext(filename)
    
    newname = indexdict[basename] + extension
    # if argument is to copy the files and then replace the filenames    
    if args.copy:
        if not os.path.exists(args.dir):
            os.makedirs(args.dir)
            
        if not os.path.exists(args.dir + "/" + newname):
            shutil.copy(fullname, os.path.join(args.dir, newname))
            print("File copied successfully.")
            print("Replacing name:" + filename + " with " + newname)
    else:
        print("Replacing name directly:" + filename + " with " + newname)
        os.rename(fullname, os.path.join(root, newname))
  


"""
rename filenames with csv file
:param args
:writes  renamed file according to csv file
:creates directory (if --copy)
"""   
def rename_file(args):
    '''
    read args.file = csv file with two columns (old_name; new_name)
    create a dictionary with the csv file: key = old_name, value = new_name
    looping through each file, and check if basename matches key of dictionary
    if matched:
      replace name with dictionary value or
      --copy: to --dir and replace name with dictionary value
    '''
    
    indexfile = pd.read_csv(args.file)
    indexdict = pd.Series(indexfile.new_name.values,index=indexfile.old_name).to_dict()
    
    # walk through each file in input directory
    for root, dirs, files in os.walk(args.input):
        for filename in files:
            fullname = os.path.join(root, filename)
            
            if args.ext:
                basename = filename.split("." + args.ext,)[0]
            else:
                basename, extension = os.path.splitext(filename)
                
            # if the basename matches the dictionary key
            if basename in indexdict.keys():
                if args.ext:
                    if fullname.endswith(args.ext):
                        rename_func(args, indexdict, fullname, filename)
                    else:
                        pass
                else:
                    rename_func(args, indexdict, fullname, filename)
                    
   
                
def main():
    parser = argparse.ArgumentParser(description="Replace multiple file names in a directory using reference csv file")
    parser.add_argument('-i','--input',help='\t\tInput directory of files to replace')
    parser.add_argument('-f','--file',help='\t\tcsv reference file <old_name>,<new_name>')
    parser.add_argument('-e','--ext',required=False, help='\t\tfilename extension to match and extract')
    parser.add_argument('-c','--copy',default=False, action="store_true", help='\t\tto copy and replace files')
    parser.add_argument('-d','--dir',required=False, help='\t\tnew directory to copy and replace files.')
  
    args = parser.parse_args()
    
    if args.copy:
        if args.dir is None:
            print("--dir required for copying files")
            sys.exit()
    
    rename_file(args)
    
    print("All Done")
    

if __name__ == "__main__":
    main()
