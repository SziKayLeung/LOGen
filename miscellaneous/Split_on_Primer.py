#!/usr/bin/env python

# taken https://github.com/Y-Lammers/Split_on_Primer/blob/master/src/Split_on_Primer.py and converted from python 2 to python 3 (indentation error)

# Split a fasta or fastq file on a set of primer sequences.
# This script can be used to seperate multiplexed sequence data
# into seperate files.

# The script requires both an input sequence file (fasta or fastq)
# and a comma seperated primer file, containing the primer
# name and primer sequence.

# Arguments: Split_on_Primer.py -s [sequence file] -p [primer file] -r -m [mismatched nucleotides allowed] -c [number of threads]

# author: Youri Lammers
# contact: youri.lammers@naturalis.nl / youri.lammers@gmail.com

# import the modules used by the script
import os, argparse, itertools, sys, multiprocessing, csv

# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description = 'Split a sequence file based on a list of primers.')

parser.add_argument('-f', '--sequence_file', metavar='Sequence file', dest='sequence', type=str,
			help='The sequence file in either fastq or fasta format.')
parser.add_argument('-p', '--primers', metavar='Primer file', dest='primer', type=str,
			help='Separated value file containing the primers. Format = primer_name,primer_sequence')
parser.add_argument('-m', '--mis', metavar='Mismatches allowed', dest='mis', type=int,
			help='The maximum number of mismatches allowed between the primer and reads (default = 0)', default=0)
parser.add_argument('-s', '--shift', metavar='Nucleotide shift allowed', dest='shift', type=int,
			help='The maximum sequence shift allowed between the primer and reads (default = 0)', default=0)
parser.add_argument('-t', '--trim', action='store_true', dest='trim',
			help='Trim the primers from the sequences after splitting.')
parser.add_argument('-d', '--delimiter', metavar='CSV delimiter', dest='delimiter', type=str,
			help='CSV delimiter used in the primer file (default = \',\')', default=',')
parser.add_argument('-c', '--cores', metavar='Number of Cores', dest='cores', type=int,
			help='The number of CPU cores the script will use (default = max number of CPUs available)', default=multiprocessing.cpu_count())
parser.add_argument('--chunk', metavar='Chunk size', dest='chunk_size', type=int,
			help='The maximum number of reads that will be loaded into the memory.\nA higher value will be faster but will take up more RAM space. (default = 25.000 * number of CPUs)', default=(25000*multiprocessing.cpu_count()))

args = parser.parse_args()

def extract_sequences():

    # open the sequence file submitted by the user, get 
    # the file format and rewind the file
    sequence_file = open(args.sequence)
    file_format = sequence_file.readline()[0]
    sequence_file.seek(0)

    # create a iterative index of all the headers
    lines = (x[1] for x in itertools.groupby(sequence_file, key=lambda line: line[0] == file_format))

    # walk through the header and obtain the sequences (and quality score if applicable)
    for headers in lines:
        header = headers.next().strip()
        if file_format == '>': sequence = [''.join(line.strip() for line in lines.next()).upper()]
        else:
            temporary_list, sequence, quality = [line.strip() for line in lines.next()], [], []

            # get the multi line sequences and break at the sequence - quality
            # seperator symbol (+)
            while len(temporary_list) > 0:
                line = temporary_list.pop(0)
                if line[0] == '+':
                    break
                sequence.append(line)
            quality = temporary_list

            # if the length of the sequences differs from the length of the
            # quality scores (because the quality line starts with a '@') get
            # the next quality line and append it to the previous one
            while len(quality) < len(sequence):
                if len(quality) == 0 and len(sequence) == 1:
                    quality.append(headers.next().strip())
                else:
                    quality += [line.strip() for line in lines.next()]
            
            # join the sequence lines and quality lines together
            sequence = [''.join(sequence).upper(), ''.join(quality)]
                
        # yield the header + sequence
        yield [header, sequence]



def extract_primers():

    # Split current primer list to a list and dictionary
    # List to keep track of # of hit per primer
    # Dictionary for easy access to output files

    # create the primer list, the list format is:
    # [primer_name, primer_squence+shift, primer_file, original_length]
    primer_list, file_dictionary = {}, {}

    # set the output dictionary (same folders as the sequence file)
    # and get the extention for the ouput files (sames as input file)
    #directory = os.path.dirname(os.path.realpath(args.sequence)) + '/'
    #directory = os.path.splitext(args.sequence)[0] + '-'
    directory = os.path.dirname(os.path.realpath(args.sequence)) + '/'
    base = os.path.splitext(os.path.basename(args.sequence))[0]
    extension = os.path.splitext(args.sequence)[1]
    directory = directory + base.replace('_', '-') + '-'

    # sanatize possible tab delimiters
    if args.delimiter == 'tab': args.delimiter = '\t'

    # walk through the primers in the primer file
    primer_file = csv.reader(open(args.primer), delimiter=args.delimiter.decode('string_escape'))
    for line in primer_file:

        # skipp comment lines in primer file
        if line[0][0] == '#': continue

        # sanatize the primer names
        line[0] = ''.join([l for l in line[0] if l.isalnum() or l in '-_'])

        # open the output_file
        output_file = open(directory + line[0] + extension, 'w')
        if line[0] not in file_dictionary:
            print('{0}{1}{2}'.format(directory, line[0], extension))

        # if the --shift argument > 0, create different primer
        # variants with the sequence shifts needed
        length = len(line[1])
        for i in range(0, args.shift + 1):
            primer_list.append([line[0], line[1][i:].upper(), length])

        # add output file to the file_dictionary
        file_dictionary[line[0]] = output_file

    # create the unsorted file
    file_dictionary['unsorted'] = open(directory + 'unsorted' + extension, 'w')
    print('{0}{1}{2}'.format(directory, line[0], extension))

    # return the primer dictionary
    return [primer_list, file_dictionary]



def compare_nuc(p, n):

	# compare the primer and sequence nucleotides
	# accounting for ambigious nucleotides in the pimer sequence

	# ambigiuous dictionary
	ambi = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
		'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
		'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']}

	# if the primer nucleotide is an ambigious character
	# check all the possible translations, return FALSE if present
	if p in ambi:
		temp = []
		for sub in ambi[p]:
			temp.append(sub != n)
		return min(temp)
	# if both nucleotides are normal, make just return the compbare boolean
	else:
		return p != n


def hamming_distance(sequence, primer):

	# calculate the hamming distance between the primer and read sequence
	# return the calculated distance
	return sum([compare_nuc(p_nuc,s_nuc) for s_nuc, p_nuc in zip(sequence, primer)])


def levenshtein_distance(sequence, primer):
	
	# calculate and return the levenshtein distance for the two sequences.
	# the levenshtein distance is only calculated if the hamming distance
	# was larger than the maximum mismatch and the --mis argument is
	# greater than 0 (default)

	previous = xrange(len(sequence) + 1)
	for pos_seq, nuc_seq in enumerate(sequence):
		current = [pos_seq + 1]
		for pos_prim, prim_seq in enumerate(primer):
			insert, delete, change = previous[pos_prim + 1]+1, current[pos_prim]+1, previous[pos_prim] + compare_nuc(prim_seq,nuc_seq)
			current.append(min(insert, delete, change))
		previous = current

	# return the distance
	return previous[-1]


def trim_primer(sequence, length):

	# trim the primer from the sequence (either in fasta or fastq format)
	# function is only used if the --trim arguments is provided
	sequence[1][0] = sequence[1][0][length:]
	if len(sequence[1]) > 1:
		sequence[1][1] = sequence[1][1][length:]
	
	# return the trimmed sequence
	return sequence


def compare_sequences(sequence, primer_list, read_shift, method, distance_results, max_mis):

	# compare the sequence to the primers

	# parse through the primers in the primer dictionary
	for primer in primer_list:

		# get the primer information
		primer_name, primer_sequence, primer_length = primer
		primer_shift = primer_length - len(primer_sequence)

		# skip same shift comparisons between the read and primers unless the shift
		# equals zero, futhermore skip shift that have been carried out before for the
		# primer - read combination. (ie read shift 2, primer shift 1 equals 
		# read shift 1, primer shift 0)
		if read_shift != 0 and primer_shift != 0: continue

		# calculate the distance with either the hamming or levenshtein method		
		if method == 'hamming': 
			distance = hamming_distance(sequence[:len(primer_sequence)], primer_sequence) + abs(read_shift - primer_shift)
		else:
			distance = levenshtein_distance(sequence[:len(primer_sequence)], primer_sequence) + abs(read_shift - primer_shift)
		
		# append the distance results if they are lower than the mis threshold
		if distance <= max_mis: 
			distance_results.append([distance, primer_name, primer_length])
			break
		
	# return the calculated distance_results list (sorted)
	return distance_results


def find_best_primer(read, primer_list, size):

	# worker thread for the distance calculations
	# this function will calculate the distance
	# between the pimers and sequence with either the
	# hamming or levenhstein method

	# list with the distance results
	distance_results = []

	# if the sequence is too short, return the read as unsorted
	if len(read[1][0]) <= size:
		return (read, 'unsorted')

	# create the sequence for each potential read shift indicated by
	# the --shift arugment
	for read_shift in range(0,args.shift+1):

		# set the shifted sequence
		sequence = read[1][0][read_shift:]

		# obtain the distance results for all primers using
		# the hamming distance method
		distance_results = compare_sequences(sequence, primer_list, read_shift, 'hamming', distance_results, args.mis)

		# if a result is obtained, jump out of the loop and return the primer
		if len(distance_results) > 0:
			break

		# If no distance could be obtained with the hamming distance method, the Levenshtein
		# method is used. Method is only used if the maximum number of mismatches is larger than 1.
		distance_results = compare_sequences(sequence, primer_list, read_shift, 'levenshtein', distance_results, args.mis)

		# again jump out if a hit is found
		if len(distance_results) > 0:
			break

	
	# pick the best match if there are multiple valid
	# primers found
	if len(distance_results) >= 1:

		# sort the distance results list (ascending based on distance)
		distance_results.sort()

		# pick the best primer
		primer = distance_results[0]

		# trim the read if --trim is selected
		if args.trim == True: read = trim_primer(read, primer[2])

		# return the read and primer
		return (read, primer[1])
	else:
		# no primer matches with the read, the read unsorted
		# will be written to the unsorted file
		return (read, 'unsorted')


def write_read(read, output_file):

	# write the read to the output_file in either fasta or fastq
	# format depending on the read
	if len(read[1]) > 1:
		output_file.write('{0}\n{1}\n+\n{2}\n'.format(read[0], read[1][0], read[1][1]))
	else:
		output_file.write('{0}\n{1}\n'.format(read[0], '\n'.join([read[1][0][i:i+60] for i in range(0, len(read[1][0]), 60)])))


def main ():

	# obtain the list with the primer names
	# sequences and mismatch shifts
	primer_list, file_dictionary = extract_primers()

	# get the largest value in the primer list + max mismatch
	size = sorted(((value[2] + args.mis) for value in primer_list), reverse=True)[0]

	# set the chunk size for the generators and worker threads
	if args.chunk_size < (args.cores*10): sub_chunk = 1
	else: sub_chunk = int(args.chunk_size/(args.cores*10))

	read_generator = extract_sequences()
	while True:	

		# create a sub iterator with a N number of reads to conserve memory
		# N is provided by the --chunk argument
		group = [list([key,chunk]) for key, chunk in itertools.islice(read_generator, args.chunk_size)]
		if group:
			
			# start a number of processes that equal the number of cores provided by --cores
			# obtain the results for each process and send them to the write function
			pool = multiprocessing.Pool(processes=args.cores)
			it = pool.imap_unordered(find_best_primer, [(read,  primer_list, size) for read in group], chunksize=sub_chunk)
			for result in it:
				write_read(result[0], file_dictionary[result[1]])
			pool.close()
			pool.join()
		else:
			break
	

if __name__ == '__main__':
	main()

