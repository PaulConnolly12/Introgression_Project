# Purpose: find the windows with the largest difference in divergence between two populations,
# and then run a finer scale divergence calculation across each window


# Run with the following command:
##python divergence_outlier_analysis.py --pop_list <pop abbreviations separated by spaces> --pop_indexes <pop indexes separated by spaces> --window_divergence_file <filename> --position_divergence_file <file> --window_starts_file <file> --window_size <size> --window_type <type> --out <path and prefix>
##python divergence_outlier_analysis.py --pop_list CO ZI --pop_index_list 1 21 --window_divergence_file test_run2_total_10.div_window --position_divergence_file test_run2_total_10.div_positions --window_starts_file test_run2_total_10.win_starts --window_size 3 --window_type filtered --out outlier_testrun_1

# Import libraries 
import argparse
import linecache
import datetime
import math
### Define functions

# Parse the window divergence values
def parse_win_div(input_window_file,pops_of_interest):
	in_window_div = []
	with open(input_window_file) as win_file:
		# Read the window file one line at a time
		for (lnum, line) in enumerate(win_file):
			div_vals = line.strip().split('\t')
			div_vals_of_interest = []
			for i in pops_of_interest:
				div_vals_of_interest += [float(div_vals[i])]
			in_window_div += [div_vals_of_interest]
	return in_window_div

# Parse the window starts
def parse_win_starts(input_start_file):
	in_window_starts = []
	with open(input_start_file) as start_file:
		# Read the window file, it's only one line long :l
		line = start_file.readline()
		starts = line.strip().split('\t')
		for start in starts:
			in_window_starts += [int(start)]
	return in_window_starts

# A temporary function to calculate the chromosome length from the inputs
def get_chrom_length(div_pos_file):
    with open(div_pos_file) as f:
    	i = -1
        for i, l in enumerate(f):
            pass
    return i + 1

# Determine regions to search through at a finer scale
def find_outliers(in_window_div,in_window_starts,outlier_proportion,chrom_len):
	num_windows = len(in_window_div)
	window_ratios = []
	for i in range(num_windows):
		ratio = in_window_div[i][0]/in_window_div[i][1]
		if math.isnan(in_window_div[i][0]):
			ratio = -1
		window_ratios.append(ratio)
	#Sort window divergence values from high to low
	sorted_window_div = sorted(window_ratios, reverse=True)
	#print(sorted_window_div)
	#Calculate how many regions to pull
	pulled_region_count = int(num_windows * outlier_proportion)
	pulled_region_count = max(pulled_region_count,1)
	#Add the outliers to a new list
	div_outliers = sorted_window_div[0:pulled_region_count]
	div_indexes = []
	#Append index values from in_window_starts to assosiated values in div_outliers
	for value in div_outliers:
		div_indexes.append(window_ratios.index(value))
	#Create lists for start and end positions of each region
	outlier_start_positions = []
	outlier_end_positions = []
	# input_window_size = in_window_starts[1] - in_window_starts[0]
	outlier_pop1_div = []
	outlier_pop2_div = []

	for index in div_indexes:
		outlier_start_positions.append(in_window_starts[index])
		if index == num_windows - 1:
			outlier_end_positions.append(chrom_len - 1)
		else:
			outlier_end_positions.append(in_window_starts[index+1] - 1)
		outlier_pop1_div.append(in_window_div[index][0])
		outlier_pop2_div.append(in_window_div[index][1])
	
	regions = (outlier_start_positions,outlier_end_positions,div_outliers,outlier_pop1_div,outlier_pop2_div)	
	
	return regions


# Pull out the masked per-position divergence values for the pops of interest
def get_reg_missing_masked(regions,input_pos_div_file,pop_indexes):
	reg_starts = regions[0]
	reg_ends = regions[1]
	#print(regions)
	reg_missing_masked_sets = []
	for (i, start) in enumerate(reg_starts):
		reg_missing_masked = [[],[]]
		for j in range(start,reg_ends[i]):
			#print(linecache.getline(input_pos_div_file,j+1))
			pos_div_line = linecache.getline(input_pos_div_file,j+1).strip().split('\t')
			value_pair = []
			#print(pos_div_line)
			for k in pop_indexes:
				if pos_div_line[k] == "N":
					value_pair.append("N")
				else:
					value_pair.append(float(pos_div_line[k]))
			#print(value_pair)
			reg_missing_masked[0].append(value_pair[0])
			reg_missing_masked[1].append(value_pair[1])
		reg_missing_masked_sets.append(reg_missing_masked)

	return reg_missing_masked_sets


# Calculate filtered windows across each specified region
def calc_filt_wind_per_region(regions,win_size,reg_missing_masked_sets):

	#print(reg_missing_masked_sets)

	region_window_div_list = []
	region_window_starts = []


	
	for i in range(len(reg_missing_masked_sets)):
		pop1_divergence = reg_missing_masked_sets[i][0]
		pop2_divergence = reg_missing_masked_sets[i][1]

		region_start = regions[0][i]
		window_starts = [region_start]
		current_win_length = 0
		pop1_filtered = []
		pop2_filtered = []

		for i in range(len(pop1_divergence)):
			# Record the current window start index if a window size of usable sites has passed
			if current_win_length == win_size:
				window_starts.append(region_start + i)
				current_win_length = 0 

			# Only transfer the divergence values for non-masked sites
			if pop1_divergence[i] != "N":
				pop1_filtered.append(pop1_divergence[i])
				pop2_filtered.append(pop2_divergence[i])
				current_win_length += 1

		region_window_starts.append(window_starts)

		window_div_list = []

		# For each population list in filtered list of counts
		for pop in [pop1_filtered,pop2_filtered]:
			# Set a start index for window
			start = 0
			# Set an end index for window that is equal to specified window size
			end = win_size
			# List that stores per-site divergence values for each window
			window_div_pop = list()

			# Perform per-site divergence calculations across each window
			for window in range(0, (((len(pop)-1)//win_size)+1)):
				# Sum all of the divergent counts across window
				sum_div = sum(pop[start:end])
				# Get per-site divergence by dividing these by the number of sites
				div_ps = sum_div/len(pop[start:end])
				# Increment start and end indices by the specified size of the window
				start += win_size
				end += win_size
				# Add this window's per-site divergence to the list
				window_div_pop.append(div_ps)

			# Add the list of window values for this population to the master list
			window_div_list.append(window_div_pop)
		# Add the list of population specific window values for this region to the master list
		region_window_div_list.append(window_div_list)

	#status = "Number of windows (window size based on called sites): " + str(len(region_window_div_list[0]))
	status = "Finished fine window divergence calcutlations per region"
	return region_window_div_list,region_window_starts,status


#rather than a list of pop with their div, just two named po
#extra for loop to go through each region
#run through divergence calcs
#get output into an output file

# Calculate fixed genomic length windows across each specified region
def calc_total_wind_per_region(regions,win_size,reg_missing_masked_sets):

	return


# Create output table composed of columns of region indexes, population indexes, window starts, window ends,
#   regions starts, region ends, and associated divergence values
def create_table(region_window_div,region_window_starts,regions,pop_names,out):

	num_entries = 0

	# Write outfile using path and name specified in input argument
	with open((out), "w") as outfile:
		headerline = "region_index\twindow_index\t"+str(pop_names[0])+"_divergence\t"+str(pop_names[1])+\
			"_divergence\t"+str(pop_names[0])+"_to_"+str(pop_names[1])+"_ratio\twindow_starts\twindow_ends\tregion_starts\tregion_ends\tregion_"\
			+str(pop_names[0])+"_divergence\tregion_"+str(pop_names[1])+"_divergence\tregion_"+str(pop_names[0])+"_to_"+str(pop_names[1])+"_ratio\n"
		outfile.write(headerline)
		for i in range(0, len(region_window_div)):
			region_start = regions[0][i]
			region_end = regions[1][i]
			region_div_ratio = regions[2][i]
			region_pop1_div = regions[3][i]
			region_pop2_div = regions[4][i]
				
			for k in range(0,len(region_window_div[i][0])):
				# Record the region/window indixes
				dataline = str(i) + "\t" + str(k) + "\t"
				# Record the window divergence data
				pop1_win_div = region_window_div[i][0][k]
				pop2_win_div = region_window_div[i][1][k]
				win_ratio = pop1_win_div/pop2_win_div
				dataline += str(pop1_win_div) +"\t" + str(pop2_win_div) + "\t" + str(win_ratio) + "\t"
				# Record the region and window start and end positions
				win_start = region_window_starts[i][k]
				win_end = 0
				if k == len(region_window_div[i][0])-1:
					win_end = region_end
				else:
					win_end = region_window_starts[i][k+1] - 1
				dataline += str(win_start) +"\t" + str(win_end) + "\t" + str(region_start) + "\t" + str(region_end) + "\t"
				# Record the region divergences
				dataline += str(region_pop1_div) +"\t" + str(region_pop2_div) + "\t" + str(region_div_ratio) + "\n"
				outfile.write(dataline)
				# Account the total number of windows
				num_entries += 1
	
	status = "Divergence data output with " + str(num_entries) + " total new windows over " + str(len(region_window_div)) + " outlier regions"
	return status


# Input arguments
def parse_args():

	# Options for input files
	parser = argparse.ArgumentParser()

	parser.add_argument('--pop_list', nargs='+', required=True,
						help='Provide the list of populations of interest')
	parser.add_argument('--pop_index_list', nargs='+', required=True,
						help='Provide the list of indexes of populations of interest')
	parser.add_argument('--window_divergence_file', required=True,
						help='Provide the window divergence file to be used in analysis')
	parser.add_argument('--position_divergence_file', required=True,
						help='Provide the per-position divergence file to be used in analysis')
	parser.add_argument('--window_starts_file', required=True,
						help='Provide the window starts file to be used in analysis')
	# parser.add_argument('--input_window_size', required=False,
	# 					help='Size of windows to calculate divergence across',
	# 					default=100000)
	parser.add_argument('--d_mel_count_files', nargs='+', required=False,
						help='Provide the list of D. mel population count files to be used in analysis')
	parser.add_argument('--d_sim_ref_count_file', required=False,
						help='Provide the name of the D. sim reference count file')
	parser.add_argument('--d_mel_count_file_path', required=False,
						help='Provide the path to these count files')
	parser.add_argument('--d_sim_count_file_path', required=False,
						help='Provide the path to this count file')
	parser.add_argument('--window_size', required=True,
						help='Size of windows to calculate divergence across')
	parser.add_argument('--window_type', required=False,
						help='Can be filtered or total depending on whether to divide window based on the number of filtered sites or the total number of sites',
						default='filtered')
	parser.add_argument('--outlier_proportion', required=False,
						help='What proportion of outliers to call and analyze per Chromosome',
						default=0.025)
	parser.add_argument('--out', required=True,
						help='Provide path and prefix of output files')

	args = parser.parse_args()

	return args



# Main function
def main():

	# Parse arguments
	args = parse_args()

	# Generate file name variables from input arguments
	#d_mel = [args.d_mel_count_file_path + f for f in args.d_mel_count_files]
	#d_sim = str(args.d_sim_count_file_path + args.d_sim_ref_count_file)
	input_window_file = str(args.window_divergence_file)
	input_pos_div_file = str(args.position_divergence_file)
	input_start_file = str(args.window_starts_file)
	out = str(args.out) + "_" + str(args.window_type) + "_" + str(args.window_size) + ".fine_div_window"
	# start = str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".fine_win_starts"
	stat = open((str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".status"), 'w')

	# Generate further input parameters
	# input_window_size = int(args.input_window_size)

	# Generate the analysis parameters
	win_type = str(args.window_type)
	size = int(args.window_size)
	outlier_proportion = float(args.outlier_proportion)
	pop_names = args.pop_list
	pop_input_indexes = [int(i) for i in args.pop_index_list]



	stat.write("Running a fine-scale divergence analysis at " + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + "\n")
	stat.write("Window size: " + str(size) + "\n" + "Window type: " + str(win_type) + "\n")

	# Parse the inputs into a list of 
	# [[pop1_win_div1,pop2_win_div1],[pop1_win_div1,pop2_win_div1],..] pre-calculated widow divergences
	in_window_div = parse_win_div(input_window_file,pop_input_indexes)
	# and [start1,start2,..] associated window starts
	in_window_starts = parse_win_starts(input_start_file)

	# print(in_window_div[0:5])
	# print(in_window_starts[0:5])

	# Identify all outlier regions
	# Add , stat line?
	chrom_len = get_chrom_length(input_pos_div_file)
	regions = find_outliers(in_window_div,in_window_starts,outlier_proportion,chrom_len)
	# stat.write(stat_line + "\n")

	# Retrieve the pre-calculated per-position divergence values
	reg_missing_masked_sets = get_reg_missing_masked(regions,input_pos_div_file,pop_input_indexes)
	# stat.write(stat_line + "\n")

	fine_div_win_list=[]
	fine_win_starts=[]
	# Determine whether total number of sites or filtered number of sites will be used to set window size
	if (win_type == 'filtered'):
		# Calculate per-site divergence within windows within regions, gives list of lists of per-window divergence in each region
		fine_div_win_list,fine_win_starts,stat_line = calc_filt_wind_per_region(regions,size,reg_missing_masked_sets)
		stat.write(stat_line + "\n")
	if (win_type == 'total'):
		fine_div_win_list,fine_win_starts,stat_line = calc_total_wind_per_region(regions,size,reg_missing_masked_sets)
		stat.write(stat_line + "\n")


	stat_line = create_table(fine_div_win_list,fine_win_starts,regions,pop_names,out)
	stat.write(stat_line + "\n")
	stat.close()
	return
		

if __name__ == '__main__':
	main()
