import os
import sys
import argparse
import logging

sys.path.append('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/0_General_py')
import read_write_fasta as rw


def log_file(output_dir, name):
    logFile = output_dir + "/" + name + '_polyApolyT.log'
    logging.basicConfig(filename = logFile,filemode = 'w',level = logging.DEBUG,
                        format = '%(asctime)s - %(levelname)s: %(message)s',\
                        datefmt = '%m/%d/%Y %I:%M:%S %p' )


def subset_poly_A_T(all_sequence_dict):
    polyT_seq = {}
    polyA_seq = {}
    nopolyA_T_seq = {}

    for trans_id, seq in all_sequence_dict.items():
        if seq.startswith("TTTTTTT"):
            polyT_seq[trans_id] = seq
        elif seq.endswith("AAAAAAA"):
            polyA_seq[trans_id] = seq
        else:
            nopolyA_T_seq[trans_id] = seq

    return(polyT_seq, polyA_seq, nopolyA_T_seq)


def main():
    parser = argparse.ArgumentParser(description="Extracting polyA and polyT sequences")
    parser.add_argument('--fa', "--input_fasta", help='\t\tFasta sequence.')
    parser.add_argument('--o_name',"--output_name", help='\t\t Output name.')
    parser.add_argument('--o_dir', "--output_dir", help='\t\tOutput path and name for list of ONT retained transcript IDs')

    args = parser.parse_args()
    log_file(args.o_dir, args.o_name)

    logging.debug("Reading fasta")
    all_sequence_dict = rw.read_inputfasta(args.fa)

    logging.debug("Subset by polyA and polyT")
    polyT_seq, polyA_seq, nopolyA_T_seq = subset_poly_A_T(all_sequence_dict)
    logging.info("Total number of reads:" + str(len(all_sequence_dict)))
    logging.info("Number of reads with polyA:" + str(len(polyA_seq)))
    logging.info("Number of reads with polyT:" + str(len(polyT_seq)))
    logging.warning("Number of reads with no polyA or polyT:" + str(len(nopolyA_T_seq)))

    rw.prepare_output(polyT_seq, args.o_dir, args.o_name + "_PolyT")
    rw.prepare_output(polyA_seq, args.o_dir, args.o_name + "_PolyA")
    rw.prepare_output(nopolyA_T_seq, args.o_dir, args.o_name + "_NoPolyA_T")
    logging.info("All Done")

if __name__ == "__main__":
    main()