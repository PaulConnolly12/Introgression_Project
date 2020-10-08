# Purpose: convert D. yakuba FASTA files into count file format

# Notes: 
## Sites with lower case or 'n' values are assumed to be missing
## This script only works for reference genome data and should not be used for population data

# Run with the following command:
## python process_counts.py --d_yak_fasta <list of yakuba FASTA files>
## python process_counts.py --d_yak_fasta yakuba_ChrX.fa yakuba_Chr2R.fa yakuba_Chr2L.fa yakuba_Chr3R.fa yakuba_Chr3L.fa

# Import libraries
import argparse
import re

# Define functions

# Input arguments
def parse_args():

	# Define script arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--d_yak_fasta', nargs='+', required=True,
						help='Provide the list of D. yakuba fasta files to be converted to count files')
	args = parser.parse_args()
	return args

# Convert FASTA to count and write to output file
def convert_to_count(file):

	# Break apart input file name
	f_name = file.replace('.', '_')
	in_name = f_name.split('_')
	# Define output file name
	out_name = 'AllCounts_yak_' + str(in_name[1]) + '.txt'
	# Open output file
	outfile = open(out_name, 'w')

	with open(file) as fasta:

		for line in fasta:
			# Skip FASTA header
			if ('>' not in line):
				# Break sequence into individual letters
				elements = list(line)

				# Write output line based on nucleotide in FASTA
				for letter in elements:
					line_to_print = ''
					if (str(letter) == 'A'):
						line_to_print = '1' + '\t' + '0' + '\t' + '0' + '\t' + '0'
					elif(str(letter) == 'C'):
						line_to_print = '0' + '\t' + '1' + '\t' + '0' + '\t' + '0'
					elif(str(letter) == 'G'):
						line_to_print = '0' + '\t' + '0' + '\t' + '1' + '\t' + '0'
					elif(str(letter) == 'T'):
						line_to_print = '0' + '\t' + '0' + '\t' + '0' + '\t' + '1'
					else:
						line_to_print = '0' + '\t' + '0' + '\t' + '0' + '\t' + '0'
					outfile.write(line_to_print + '\n')
	outfile.close()

# Main function
def main():

	# Parse inputs
	args = parse_args()
	# List of input files
	files = args.d_yak_fasta

	# Write output for each input file
	for f in files:
		convert_to_count(file=f)

# Run main function
if __name__ == '__main__':
	main()

				