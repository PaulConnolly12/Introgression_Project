# Purpose: calculate per-site divergence between D. sim for arbitrary number of D. mel populations

# Notes: 
## Only include sites with at least one called individual for every population in div calculation
## Need to keep track of how many sites are being used in total
## [TO DO] output the genomic coordinates for windows
## [TO DO] include the other parameter values used in fnames output

# Run with the following command:
## python divergence.py --d_mel_count_files <list of files> --d_sim_ref_count_file <file> --d_mel_count_file_path <path to files> --d_sim_count_file_path <path to file> --window_size <int> --window_type <filtered or total> --out <path + filename>
## python divergence.py --d_mel_count_files AllCounts_B_Chr2R.txt AllCounts_CO_Chr2R.txt AllCounts_EA_Chr2R.txt AllCounts_EF_Chr2R.txt AllCounts_EG_Chr2R.txt AllCounts_FR_Chr2R.txt AllCounts_GA_Chr2R.txt AllCounts_GU_Chr2R.txt AllCounts_I_Chr2R.txt AllCounts_KF_Chr2R.txt AllCounts_Kenya_Chr2R.txt AllCounts_MW_Chr2R.txt AllCounts_NG_Chr2R.txt AllCounts_N_Chr2R.txt AllCounts_RAL_Chr2R.txt AllCounts_RG_Chr2R.txt AllCounts_SD_Chr2R.txt AllCounts_SP_Chr2R.txt AllCounts_T_Chr2R.txt AllCounts_WAF_Chr2R.txt AllCounts_W_Chr2R.txt AllCounts_ZI_Chr2R.txt --d_sim_ref_count_file AllCounts_sim_Chr2R.txt --d_mel_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/ --d_sim_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/ --window_size 10000 --window_type filtered --out test_filtered
## python divergence.py --d_mel_count_files toydata_CO_Chr2R.txt toydata_FR_Chr2R.txt --d_sim_ref_count_file toydata_sim_Chr2R.txt --d_sim_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/ --d_mel_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/ --window_size 10 --window_type total --out /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_run 

# Import libraries 
import argparse
import linecache

# Define functions

# Input arguments
def parse_args():

	# Options for input files
	parser = argparse.ArgumentParser()

	parser.add_argument('--d_mel_count_files', nargs='+', required=True,
						help='Provide the list of D. mel population count files to be used in analysis')
	parser.add_argument('--d_sim_ref_count_file', required=True,
						help='Provide the name of the D. sim reference count file')
	parser.add_argument('--d_mel_count_file_path', required=True,
						help='Provide the path to these count files')
	parser.add_argument('--d_sim_count_file_path', required=True,
						help='Provide the path to this count file')		
	parser.add_argument('--window_size', required=True,
						help='Size of windows to calculate divergence across')
	parser.add_argument('--window_type', required=True,
						help='Can be filtered or total depending on whether to divide window based on the number of filtered sites or the total number of sites')
	parser.add_argument('--out', required=True,
						help='Provide path and prefix of output files')

	args = parser.parse_args()
	return args


# Calculate proportion of melanogaster genotypes that don't match simulans reference
def div_calc(d_mel, d_sim):
	#print(d_sim[0:10])
	#print(d_mel[0:10])

	# Define a list of lists
	pop_list = list()

	# Iterate through the list of d_mel population count files
	for f in d_mel:
		# Create list where each element is divergence at that site index
		div_list = list()
		
		# Open each melanogaster population file
		with open(f) as mel_file:
			# Read d_mel file one line at a time
			for (lnum, line) in enumerate(mel_file):
				# Extract elements as int in matching d_sim line (note: linecache starts counting at 1)
				sim_line = list(map(int, ((linecache.getline(d_sim, (lnum + 1))).split())))
				# Extract elements as int in d_mel line for current file
				mel_line = list(map(int, (line.split())))
				# Create counter to track the number of non-ref genotypes for each site in d_mel file
				div_cnt = 0
				# print(sim_line)
				# print(mel_line)
				#print(list(sim_line))
				#print(list(mel_line))

				# Add 'N' to list at sites where either the sim ref or the mel pop are not called
				if ((sum(sim_line) == 0) or (sum(mel_line) == 0)):
					# Add a null value to the list entry for this site
					div_list.append('N')
				
				# If there are called genotypes at a site
				else:
					# Check each base-pair column
					for i in range(0, len(list(sim_line))):
						# If there are non-ref genotypes
						if ((sim_line[i] != mel_line[i]) and (mel_line[i] > 0) and (sim_line[i] == 0)):
							# Count how many non-ref genotypes there are
							div_cnt += mel_line[i]

					# To get the proportion of non-ref genotypes, divide by the total number of d_mel called genotypes
					div_prop = (float(div_cnt) / sum(mel_line))
					# Add the proportion of divergent genotypes to the list entry for this site
					div_list.append(div_prop)

		# Create a list of lists where each element is a list of population divergence values
		pop_list.append(div_list)

	status = "Calculated divergent counts for " + str(len(pop_list)) + " populations" + "\n" + "Number of sites processed for each population: " + str(len(pop_list[0]))
	return pop_list,status


# Identify sites where at least one population has no called genotypes and provide indices and the number of populations with missing data
def missing_data(pop_list):

	# Create a list of indices that contain missing data in at least one population
	miss_idx = list()
	# Create a list of the number of populations for which this site had missing data
	miss_cnt_list = list()

	# Loop through each site 
	for i in range(0, len(pop_list[0])):
		# Keep track of the number of populations that have missing data at this site
		n_cnt = 0

		# Look for any 'N's at this site across populations
		for pop in pop_list:
			# Set missing data to true if an 'N' is encountered
			if (pop[i] == 'N'):
				n_cnt += 1

		# Append this site index to list of indices if it contains missing data
		if (n_cnt>0):
			miss_idx.append(i)
			miss_cnt_list.append(n_cnt)
	status = str(len(miss_idx)) + " sites did not meet sample threshold due to missing data"
	return miss_idx, miss_cnt_list, status


# Mask sites where at least one population has no data in all population lists
def mask_missing(pop_list, miss_idx):
	
	# Create a list of lists that will store masked counts
	filtered_pop_list = list()

	# Replace sites with missing data in at least one population to 'N' in all populations
	for pop in pop_list:
		for index in miss_idx:
			pop[index] = 'N'

		filtered_pop_list.append(pop)

	status = str(len(list(filter(lambda x: x != 'N', filtered_pop_list[0])))) + " sites were retained after masking"
	return filtered_pop_list,status


# Output the per-site divergence across all (available) sites for each population
def full_div(missing_masked):
	
	# List lists where masked missing sites have been removed from count data
	filter_list = list()
	# List where each element is total per-site divergence in one population
	full_div_list = list()
	# Sites used to calculate divergence
	sites_used = list()

	# For each population list of counts
	for pop in missing_masked:
		# Remove sites with 'N'
		filter_pop = list(filter(lambda x: x != 'N', pop))
		# Create filtered data list
		filter_list.append(filter_pop)

	# For each population in filtered list of counts
	for pop in filter_list:
		# Find total per-site divergence
		full_div = sum(pop)/len(pop)
		# Create a list where each element is total per-site divergence for single population
		full_div_list.append(full_div)
		# Record number of site used to calculate divergence in this population
		sites_used.append(len(pop))

	status = "Genome-wide per-site divergence calculated for " + str(len(sites_used)) + "\n" + "Number of sites used in genome-wide per-site divergence calculation: " + str(sum(sites_used)/len(sites_used))
	return full_div_list,status


# Calculate per-site divergence across windows sizes that include only valid sites (number of sites in each window will be the same) 
def per_site_filtered(size, missing_masked):

	# List of lists where masked missing sites have been removed from count data
	filter_list = list()
	# Create a list that stores each population list of window-based divergence values
	window_div_list = list()
	# Create a list that stores the window starts
	window_starts = list()

	# # For each population list of counts, remove null sites before breaking into windows
	# for pop in missing_masked:
	# 	# Remove sites with 'N'
	# 	filter_pop = filter(lambda x: x != 'N', pop)
	# 	# Create filtered data list
	# 	filter_list.append(filter_pop)

	# Populate the new filtered list of divergence values without 'N' values for each population list
	for i in range(len(missing_masked)):
		filter_list.append([])
	# Simultaneosly populate the list of start index values for each window
	window_starts.append(0)
	current_win_length = 0
	for i in range(len(missing_masked[0])):
		# Record the current window start index if a window size of usable sites has passed
		if current_win_length == size:
			window_starts.append(i)
			current_win_length = 0
		# Only transfer the divergence values for non-masked sites
		if missing_masked[0][i] != "N":
			for j in range(len(filter_list)):
				filter_list[j].append(missing_masked[j][i])
			current_win_length += 1
	
	# For each population list in filtered list of counts
	for pop in filter_list:
		# Set a start index for window
		start = 0
		# Set an end index for window that is equal to specified window size
		end = size
		# List that stores per-site divergence values for each window
		window_div_pop = list()

		# Perform per-site divergence calculations across each window
		for window in range(0, ((len(pop)//size)+1)):
			# Sum all of the divergent counts across window
			sum_div = sum(pop[start:end])
			# Get per-site divergence by dividing these by the number of sites
			div_ps = sum_div/len(pop[start:end])
			# Increment start and end indices by the specified size of the window
			start += size
			end += size
			# Add this window's per-site divergence to the list
			window_div_pop.append(div_ps)

		# Add the list of window values for this population to the master list
		window_div_list.append(window_div_pop)

	status = "Number of windows (window size based on called sites): " + str(len(window_div_list[0]))
	return window_div_list,window_starts,status


# Calculate per-site divergence across windows sizes that include missing sites (number of sites in each window will vary) 
def per_site_total(size, missing_masked):

	# Create a list that stores each population list of window-based divergence values
	window_div_list = list()

	# List of called sites in each window in each population 
	called_sites_list = list()

	#Calculate the window start indexes
	window_starts = list()
	current_start = 0
	while current_start < len(missing_masked[0]):
		window_starts.append(current_start)
		current_start += size
	


	# For each population list in unfiltered list of counts, set window size before removing null values
	for pop in missing_masked:
		# Set start index
		start = 0
		# Set an end index for window that is equal to specified window size
		end = size
		# List that stores per-site divergence values for each window
		window_div_pop = list()

		# For each window
		for window in range(0, ((len(pop)//size)+1)):
			# Remove sites within the specified window that contain 'N'
			filter_win = list(filter(lambda x: x != 'N', pop[start:end]))

			# If the window still contains sites after filtering
			if (len(filter_win) != 0):
				# Sum the divergence counts across sites in the window
				sum_div = sum(filter_win)
				# Divide by the filtered length of the window to get per-site divergence
				div_ps = sum_div/len(filter_win)
				# Append this window divergence values to the list for this population
				window_div_pop.append(div_ps)
				# Append the number of called sites to the list for this population
				called_sites_list.append(len(filter_win))
			# If the window contains no sites after filtering, store 'NaN', Not a Number
			else:
				window_div_pop.append("NaN")
				called_sites_list.append(0)
			# Increment start and end indices by the size of the window
			start += size
			end += size

		# Append list of window-based divergences for single population to master list
		window_div_list.append(window_div_pop)

	# Calculate mean and variance of number of retained sites in windows
	mean_sites = float(sum(called_sites_list))/len(called_sites_list)
	var_sites = sum((i - mean_sites) ** 2 for i in called_sites_list)/len(called_sites_list)

	status = "Number of windows (window size based on total sites): " + str(len(window_div_list[0])) + "\n" + "Average number of called sites in each window: " + str(mean_sites) + "\n" + "Variance in number of called sites in each window: " + str(var_sites)
	return window_div_list,window_starts,status


# Create output table where column is population and row is window
def create_table(window_div, out):

	# Write outfile using path and name specified in input argument
	outfile = open((str(out)), "w")
	# Initialize list which will hold row values (same window in different populations)
	rows = list()

	# Transform the list of windows in each population to list of populations for each window
	# For every window
	for i in range(0, len(window_div[0])):
		# Create a list of that windows's divergence value in each population
		col_vals = list()

		# For each element (population) in the list of population windows
		for win in window_div:
			# Add it to list that holds the same window value in each population
			col_vals.append(win[i])

		rows.append(col_vals)
		# Print this row to the output file
		line_to_print = "\t".join([str(x) for x in col_vals]) + "\n"
		outfile.write(line_to_print)
	
	status = "Window based divergence outfile contains " + str(len(rows)) + " rows and " + str(len(rows[0])) + " columns"
	return status


# Create output table where first column is the index (counting from 0) of sites that were excluded and second column is the number of populations that did not meet threshold at that site
def create_missing(miss_idx, miss_cnt, missing):

	# Create output with the path and filename specified in input arguments
	outfile = open((str(missing)), "w")

	# Print each element in the missing indices list and missing counts list side by side with a tab in between
	for i in range(0, len(miss_idx)):
		line_to_print = "\t".join([str(miss_idx[i]), str(miss_cnt[i])]) + "\n"
		outfile.write(line_to_print)
	outfile.close()


# Create output file with an ordered list of genome-wide per-site divergence values for each population
def create_div(full_div_list, div):

	# Create output with the path and filename specified in input arguments
	outfile = open((str(div)), "w")
	# Print out the list containing genome-wide per-site divergence value for each population
	line_to_print = "\t".join([str(x) for x in full_div_list]) + "\n"
	outfile.write(line_to_print)
	outfile.close()


# Create output file with the ordered list of population files used in this analysis
def create_fname_file(filenames, file):

	# Create output with the path and filename specified in input arguments
	outfile = open((str(file)), "w")
	# Print out the list of filename arguments for the d_mel files
	line_to_print = "\t".join([str(x) for x in filenames])
	outfile.write(line_to_print)
	outfile.close()

#Create an output file with the list of window starts
def create_starts(window_starts,file_name):
	outfile = open((str(file_name)),"w")
	line_to_print = "\t".join([str(x) for x in window_starts])
	outfile.write(line_to_print)
	outfile.close()

#Creste an output file with the divergence positions
def create_divergence_position(missing_masked,file_name):
	with open((str(file_name)),"w") as outfile:
		for i in range(len(missing_masked[0])):
			line_to_print = "\t".join([str(missing_masked[j][i]) for j in range(len(missing_masked))]) + "\n"
			outfile.write(line_to_print)

# Main function
def main():

	# Parse arguments
	args = parse_args()

	# Assign input arguments to file name variables
	d_mel = [args.d_mel_count_file_path + f for f in args.d_mel_count_files]
	d_sim = str(args.d_sim_count_file_path + args.d_sim_ref_count_file)
	out = str(args.out) + "_" + str(args.window_type) + "_" + str(args.window_size) + ".div_window"
	div = str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".div_full"
	missing = str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".missing"
	start = str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".win_starts"
	stat = open((str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".status"), 'w')
	win_type = str(args.window_type)
	size = int(args.window_size)
	position = str(args.out) + "_" + str(args.window_type) + "_" + str(args.window_size) + ".div_positions"
	stat.write("Window size: " + str(size) + "\n" + "Window type: " + str(win_type) + "\n")

	# Create list of lists of divergence counts at each site in each population
	pop_list,stat_line = div_calc(d_mel=d_mel, d_sim=d_sim)
	stat.write(stat_line + "\n")

	# List of indices at which one or more population file (d_mel or d_sim) have 0 called individuals and counts of the number of populations that have missing data for each excluded site
	miss_idx,miss_cnt,stat_line = missing_data(pop_list=pop_list)
	stat.write(stat_line + "\n")

	# List of lists of divergence counts at each site in each population (where excluded sites are masked in all populations with an 'N')
	missing_masked,stat_line = mask_missing(pop_list=pop_list, miss_idx=miss_idx)
	stat.write(stat_line + "\n")

	# Determine whether total number of sites or filtered number of sites will be used to set window size
	if (win_type == 'total'):
		# Calculate per-site divergence within window, gives list of lists of per-window divergence in each population
		window_div,window_starts,stat_line = per_site_total(size=size, missing_masked=missing_masked)
		stat.write(stat_line + "\n")
	if (win_type == 'filtered'):
		window_div,window_starts,stat_line = per_site_filtered(size=size, missing_masked=missing_masked)
		stat.write(stat_line + "\n")

	# List of the genome-wide per-site divergence values for each population
	full_div_list,stat_line = full_div(missing_masked=missing_masked)
	stat.write(stat_line + "\n")

	# Create outputs
	stat_line = create_table(window_div=window_div, out=out)
	stat.write(stat_line + "\n")
	create_missing(miss_idx=miss_idx, miss_cnt=miss_cnt, missing=missing)
	create_div(full_div_list=full_div_list, div=div)
	create_fname_file(filenames=args.d_mel_count_files, file=(str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".names"))

	create_divergence_position(missing_masked=missing_masked, file_name=position)

	create_starts(window_starts=window_starts,file_name=start)

if __name__ == '__main__':
	main()
