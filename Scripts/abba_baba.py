# Purpose: compute Patterson's D statistic for a specified set of populations

# Notes:
## method 1) frequency-based: set ref frequencies to 0 or 1
## method 2) count-based: choose allele type by sampling from population frequencies 
## [TO DO] allow for window-based calculation

# Run with the following command:
## python abba_baba.py --p1 <p1 file> --p2 <p2 file> --p3 <p3 file> --pO <pO file> --window_opt <"on","off"> --size <int> --out <outfile prefix> --path_in <input path> --path_out <output_path> --window_type <"total", "retained">
## python abba_baba.py --p1 toydata_CO_Chr2R.txt --p2 toydata_FR_Chr2R.txt --p3 toydata_sim_Chr2R.txt --pO toydata_yak_Chr2R.txt --window_opt on --size 10 --out chr2_test --path_in /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/ --path_out /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/ --window_type total

# Import libraries 
import argparse
import linecache
import random

# Define functions

# Input arguments
def parse_args():

	# Options for input files
	parser = argparse.ArgumentParser()
	parser.add_argument('--p1', required=True,
						help='Population 1 in ABBA BABA test (e.g. CO)')
	parser.add_argument('--p2', required=True,
						help='Population 2 in ABBA BABA test (e.g. ZI')
	parser.add_argument('--p3', required=True,
						help='Population 3 in ABBA BABA test (e.g. D_simulans)')
	parser.add_argument('--pO', required=True,
						help='Outgroup population in ABBA BABA test (e.g. D_yakuba)')
	parser.add_argument('--window_opt', required=True,
						help='Specify "off" if using entire chromosome segment or "on" if using specific window')	
	parser.add_argument('--window_type', required=False,
						help='Specify "total" if setting window size based on total genomic sites'
						'or "retained" if setting window size based on sites retained for D calculation')	
	parser.add_argument('--size', required=False,
						help='Size of windows to compute D across')
	parser.add_argument('--out', required=True,
						help='Provide prefix for output files')
	parser.add_argument('--path_in', required=True,
						help='Provide path for input files')
	parser.add_argument('--path_out', required=True,
						help='Provide path for output files')

	args = parser.parse_args()
	return args

# Determine which allele types are at each site
def allele_types(site):

	# Create a list to hold allele types
	alleles = list()
	# Add allele A to list if there are genotypes in first col
	if (site[0] != 0):
		alleles.append('A')
	# Add allele C to list if there are genotypes in second col
	if (site[1] != 0):
		alleles.append('C')
	# Add allele G to list if there are genotypes in third col
	if (site[2] != 0):
		alleles.append('G')
	# Add allele T to list if there are genotypes in fourth col
	if (site[3] != 0):
		alleles.append('T')
	##print('alleles:\n' + str(alleles))
	return alleles

# Create a list of A, C, G, T genotype proportions in each of the four populations
def allele_list_fun(pop_list):

	# Initialize four empty lists
	a = list()
	c = list()
	g = list()
	t = list()

	# For each site in each population
	for pop in pop_list:
		# Create a list for the proportion of genotypes with A
		a.append(float(pop[0])/sum(pop))
		# Create a list for the proportion of genotypes with C
		c.append(float(pop[1])/sum(pop))
		# Create a list for the proportion of genotypes with G
		g.append(float(pop[2])/sum(pop))
		# Create a list for the proportion of genotypes with T
		t.append(float(pop[3])/sum(pop))

	# Store this genotype data in a list of lists
	allele_list = [a, c, g, t]
	##print('allele_list:\n' + str(allele_list))
	return allele_list

# Determine which two allele type lists represent the ancestral and derived states
def anc_vs_der(allele_list):

	# Initialize two empty lists to represent each population's ancestral or derived genotype proportions
	ancestral,derived = ([], ) * 2

	# For each nucleotide type
	for nuc_type in allele_list:
		# If there are genotypes for the allele type
		if (sum(nuc_type) != 0):
			# It is ancestral if contained in the outgroup
			if (nuc_type[3] == 1.0):
				ancestral = nuc_type
			# Derived otherwise
			else:
				derived = nuc_type
	##print('ancestral:\n' + str(ancestral) + '\n' + 'derived:\n' + str(derived))
	return ancestral,derived

# Perform the count-based ABBA BABA method
def count_based(ancestral,derived):

	# Initialize an empty list to hold the pattern
	pattern= list()
	# Initialize an empty list to track the number of sites whose draw outcome doesn't match ABBA or BABA
	count_rejects = list()
	# Initialize variables to hold ABBA vs BABA value (in this case 0 or 1)
	abba = 0
	baba = 0

	# Count-based method to ABBA BABA
	for i in range(0,len(ancestral)):
		# Assign A to pattern if population has ancestral allele
		if (ancestral[i] == 1.0):
			pattern.append('A')
		# Assign B to pattern if population has derived allele
		elif (derived[i] == 1.0):
			pattern.append('B')
		# Assign allele state in the event that both alleles are segregating in population
		else:
			# Pick a random float between 0 and 1.0
			draw = random.random()
			##print('draw:\n' + str(draw))
			# If derived frequency is greater than ancestral frequency
			if (derived[i] > ancestral[i]):
				# Assign state to ancestral if draw is within ancestral frequency
				if (draw < ancestral[i]):
					pattern.append('A')
				# Assign state to derived if draw is beyond ancestral frequency
				elif (draw > ancestral[i]):
					pattern.append('B')
			# If derived frequency is smaller than ancestral frequency
			elif (derived[i] < ancestral[i]):
				# Assign state to derived if draw is within derived frequency
				if (draw < derived[i]):
					pattern.append('B')
				# Assign state to ancestral if draw is beyond derived frequency
				elif (draw > derived[i]):
					pattern.append('A')
			# If the ancestral and derived allele frequencies are equal
			else:
				# Assign to ancestral if draw is within ancestral frequency
				if (draw < ancestral[i]):
					pattern.append('A')
				# Assign to derived if draw is beyond ancestral frequency
				elif (draw > ancestral[i]):
					pattern.append('B')
	
	# Turn the ABBA or BABA pattern into a word
	string = ''
	word = string.join(pattern)
	##print('word:\n' + str(word))

	# Add counts of ABBA vs BABA to a list
	if (word == 'ABBA'):
		abba = 1
	elif (word == 'BABA'):
		baba = 1
	# Keep track of sites whose draw outcome does not match ABBA or BABA pattern
	else:
		count_rejects.append('NP')
	
	##print('abba count:\n' + str(abba) + '\n' + 'baba count:\n' + str(baba))
	return abba,baba,count_rejects

# Perform frequency-based calculation
def freq_based(derived):

	# Calculate ABBA pattern weight
	abba = (1 - float(derived[0])) * float(derived[1]) * float(derived[2]) * (1 - float(derived[3]))
	# Calculate BABA pattern weight
	baba = float(derived[0]) * (1 - float(derived[1])) * float(derived[2]) * (1 - float(derived[3]))
	##print('abba freq:\n' + str(abba) + '\n' + 'baba freq:\n' + str(baba))
	return abba,baba

# Calculate the D statistic from ABBA and BABA lists
def calc_d(abba,baba):

# [TO DO] allow for window-based calculation of statistic
# This would require inserting null values into abba and baba lists for sites that aren't retained

	# Sum the ABBA vs BABA pattern values
	sum_abba = sum(abba)
	sum_baba = sum(baba)
	# Calculate the numerator of D by dividing the sum of ABBAs by the sum of BABAs
	num = sum_abba - sum_baba
	# Calculate the denominator of D by summing ABBA and BABA 
	den = sum_abba + sum_baba
	# Compute D
	d_val = float(num)/float(den)
	# Write a status line
	status = 'The sum of ABBA patterns is: ' + str(sum_abba) + '\n' + 'The sum of BABA patterns is: ' + str(sum_baba) + '\n' + 'The value of Pattersons D is: ' + str(d_val)
	return status

# Output the number of sites rejected for (1) missing data (2) not biallele (3) P3 = PO (4) doesn't match ABBA or BABA after count-based draw
def reject(total,count,num_tot,num_ret):

	# Initialize variables to hold counts of each reject type
	missing = 0
	biallelic = 0
	same_3O = 0
	same_12 = 0

	# Evaluate what the nature of each reject type is from total list
	for i in total:
		if (str(i) == 'M'):
			missing += 1
		elif (str(i) == 'B'):
			biallelic += 1
		elif (str(i) == 'S_3O'):
			same_3O += 1
		elif (str(i) == 'S_12'):
			same_12 += 1

	# Count the number of occurances for reject type 4
	count_no_match = len(count)
	# Write a status line
	status = str(num_tot) + ' sites analyzed in total\n' + str(num_ret) + ' sites retained for analysis\n' + str(missing) + ' sites were excluded due to missing data\n' + str(biallelic) + ' sites were excluded for not being biallelic\n' + str(same_3O) + ' sites were excluded because P3 genotype matched PO genotype\n' + str(same_12) + ' sites were excluded because P1 genotype matched P2 genotype\n' + str(count_no_match) + ' sites were excluded from count-based analysis because pattern did not conform to ABBA or BABA\n'
	return status

# Main function
def main():

	# Parse arguments
	args = parse_args()

	# Assign input arguments to path variables
	path_in = str(args.path_in)
	path_out = str(args.path_out)

	# Assign input arguments to file name variables
	file_names = [args.p1, args.p2, args.p3, args.pO]
	p1_file,p2_file,p3_file,pO_file = [path_in + str(f) for f in file_names]
	file_list = [p1_file, p2_file, p3_file, pO_file]
	outfile = open((path_out + str(args.out)), 'w')

	# Creating site lists for each of the four populations
	with open(p1_file) as p1:

		# Initialize ABBA and BABA lists for the frequency-based method and count-based method
		# [TO DO] for window based analysis, you would want these lists to hold null values for sites that are excluded for analysis
		abba_freq = list()
		baba_freq = list()
		abba_count = list()
		baba_count = list()

		# Initialize a list for rejected sites
		rejects = list()

		# Initialize a list for rejected sites by the count based method (draw outcome doesn't match ABBA or BABA)
		count_rejects = list()

		# Keep track of the number of sites being analyzed
		num_sites = 0

		# Keep track of the number of sites retained for analysis
		retained = 0

		# Read p1 file one line at a time
		for lnum, line in enumerate(p1):
			
			# Split the site's line into elements for p1
			p1_site = map(int, (line.split()))
			# Extract elements as int in matching p2 line (note: linecache starts counting at 1)
			p2_site = map(int, ((linecache.getline(p2_file, (lnum + 1))).split()))
			# Extract elements as int in matching p3 line
			p3_site = map(int, ((linecache.getline(p3_file, (lnum + 1))).split()))
			# Extract elements as int in matching pO line
			pO_site = map(int, ((linecache.getline(pO_file, (lnum + 1))).split()))

			# Only consider sites where all populations have data
			if ((sum(p1_site) and sum(p2_site) and sum(p3_site) and sum(pO_site)) != 0):
				
				# Make a list of populations
				pop_list = [p1_site, p2_site, p3_site, pO_site]
				# Create a list of allele types at site in each population
				p1_allele,p2_allele,p3_allele,pO_allele = [allele_types(site=p) for p in pop_list]		
				# Find the union of allele types at site
				total = list(set(p1_allele) | set(p2_allele) | set(p3_allele) | set(pO_allele))
				##print('total:\n' + str(total))

				# Only consider sites that are biallelic (for same two alleles in all populations)
				if (len(total) <= 2):

					# Only consider sites where P3 does not have same allele as PO
					if (set(pO_allele) != set(p3_allele)):

						# Only consider sites where P1 and P2 do not have the same allele
						if (set(p1_allele) != set(p2_allele)):
						
							# Initialize three empty lists
							ancestral,derived,pattern = ([], ) * 3
							# Create list of each population's A, C, G, T genotype proportions
							allele_list = allele_list_fun(pop_list=pop_list)
							# Determine which two allele type lists represent the ancestral and derived states
							ancestral,derived = anc_vs_der(allele_list=allele_list)
							
							# Perform the count-based ABBA BABA method
							abba,baba,count_rejects = count_based(ancestral=ancestral,derived=derived)
							# Append these values to the count-based ABBA and BABA lists
							abba_count.append(abba)
							baba_count.append(baba)
							
							# Perform the frequency-based ABBA BABA method
							abba,baba = freq_based(derived=derived)
							# Append these values to the frequency-based ABBA and BABA lists
							abba_freq.append(abba)
							baba_freq.append(baba)

							# Keep track of the number of sites retained
							retained += 1
						
						# Keep track of sites excluded for P3 and PO having 'S'ame allele
						else:
							rejects.append('S_12')
					# Keep track of sites excluded for P3 and PO having 'S'ame allele
					else:
						rejects.append('S_3O')
				# Keep track of sites excluded for not being 'B'iallelic
				else:
					rejects.append('B')
			# Keep track of sites excluded fo having 'M'issing data
			else:
				rejects.append('M')

			# Count the number of sites being analyzed
			num_sites += 1

	# Calculate the D statistic for both the count-based and frequency-based methods
	stat_count = calc_d(abba=abba_count,baba=baba_count)
	stat_freq = calc_d(abba=abba_freq,baba=baba_freq)
	stat_reject = reject(total=rejects,count=count_rejects,num_tot=num_sites,num_ret=retained)

	# Write the output
	outfile.write('Count-based test:' + '\n')
	outfile.write(stat_count + '\n')
	outfile.write('Frequency-based test:' + '\n')
	outfile.write(stat_freq + '\n')
	outfile.write(stat_reject)
	outfile.close()

if __name__ == '__main__':
	main()