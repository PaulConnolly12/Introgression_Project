# Purpose: find the windows with the largest difference in divergence between two populations,
# and then run a finer scale divergence calculation across each window


# Run with the following command:
## python divergence_outlier_analysis.py --d_mel_count_files <list of files> --d_sim_ref_count_file <file> --d_mel_count_file_path <path to files> --d_sim_count_file_path <path to file> --window_size <int> --window_type <filtered or total> --out <path + filename>
## python divergence_outlier_analysis.py --d_mel_count_files AllCounts_CO_Chr2R.txt AllCounts_ZI_Chr2R.txt --d_sim_ref_count_file AllCounts_sim_Chr2R.txt --d_sim_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/ --d_mel_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/ --window_size 100 --out /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/PFC/ 

# Import libraries 
import argparse
import linecache
import datetime

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



# Determine regions to search through at a finer scale
def find_outliers(in_window_div,in_window_starts,outlier_proportion):
	#Sort window divergence values from high to low
	sorted_window_div = sorted(in_window_div, reverse=True)
	#Calculate how many regions to pull
	pulled_region_count = int(len(in_window_div) * outlier_proportion)
	#Add the outliers to a new list
	div_outliers = sorted_window_div[0:pulled_region_count]
	div_indexes = []
	#Append index values from in_window_starts to assosiated values in div_outliers
	for value in div_outliers:
		div_indexes.append(in_window_div.index(value))
	#Create lists for start and end positions of each region
	outlier_start_positions = []
	outlier_end_positions = []
	input_window_size = in_window_starts[1] - in_window_starts[0]

	for index in div_indexes:
		outlier_start_positions.append(in_window_starts[index])
		outlier_end_positions.append(in_window_starts[index] + input_window_size)
	
	regions = (outlier_start_positions,outlier_end_positions)	
	
	return regions


# Pull out the masked per-position divergence values for the pops of interest
def get_reg_missing_masked(regions,input_pos_div_file,pop_indexes):
	(reg_starts,reg_ends) = regions

	reg_missing_masked_sets = []
	for (i, start) in enumerate(reg_starts):
		reg_missing_masked = [[],[]]
		for j in range(start,reg_ends[i]):
			pos_div_line = linecache.getline(input_pos_div_file,j).strip().split('\t')
			value_pair = []
			for k in pop_indexes:
				if pos_div_line[k] == "N":
					value_pair.append("N")
				else:
					value_pair.append(float(pos_div_line[k]))
			reg_missing_masked[0].append(value_pair[0])
			reg_missing_masked[1].append(value_pair[1])
		reg_missing_masked_sets.append(reg_missing_masked)

	return reg_missing_masked_sets


# Calculate filtered windows across each specified region
def calc_filt_wind_per_region(regions,win_size,reg_missing_masked_sets):


	return


# Calculate fixed genomic length windows across each specified region
def calc_total_wind_per_region(regions,win_size,reg_missing_masked_sets):

	return



# Input arguments
def parse_args():

	# Options for input files
	parser = argparse.ArgumentParser()

	parser.add_argument('--window_pop_list', nargs='+', required=False,
						help='Provide the list of populations used in generating the ')
	parser.add_argument('--window_divergence_file', required=True,
						help='Provide the window divergence file to be used in analysis')
	parser.add_argument('--position_divergence_file', required=True,
						help='Provide the per-position divergence file to be used in analysis')
	parser.add_argument('--window_starts_file', required=True,
						help='Provide the window starts file to be used in analysis')
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
	parser.add_argument('--window_type', required=False,
						help='Can be filtered or total depending on whether to divide window based on the number of filtered sites or the total number of sites',
						default='filtered')
	parser.add_argument('--out', required=True,
						help='Provide path and prefix of output files')

	args = parser.parse_args()
	return args

# Main function
def main():

	# Parse arguments
	args = parse_args()

	# Assign input arguments to file name variables
	d_mel = [args.d_mel_count_file_path + f for f in args.d_mel_count_files]
	d_sim = str(args.d_sim_count_file_path + args.d_sim_ref_count_file)
	input_window_file = str(args.window_divergence_file)
	input_pos_div_file = str(args.position_divergence_file)
	input_start_file = str(args.window_starts_file)
	out = str(args.out) + "_" + str(args.window_type) + "_" + str(args.window_size) + ".fine_div_window"
	start = str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".fine_win_starts"
	stat = open((str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".status"), 'w')
	win_type = str(args.window_type)
	size = int(args.window_size)

	stat.write("Running a fine-scale divergence analysis at " + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + "\n")
	stat.write("Window size: " + str(size) + "\n" + "Window type: " + str(win_type) + "\n")

	# Manually define the indexes of the population pair of interest
	pops_of_interest = [1,21]

	# Parse the inputs into a list of 
	# [[pop1_win_div1,pop2_win_div1],[pop1_win_div1,pop2_win_div1],..] pre-calculated widow divergences
	in_window_div = parse_win_div(input_window_file,pops_of_interest)
	# and [start1,start2,..] associated window starts
	in_window_starts = parse_win_starts(input_start_file)

	print(in_window_div[0:5])
	print(in_window_starts[0:5])

	# Manually define the proportion of the windows to pick as outliers
	outlier_proportion = 0.025

	# Identify all outlier regions
	# Add , stat line?
	regions = find_outliers(in_window_div,in_window_starts,outlier_proportion)
	# stat.write(stat_line + "\n")

	# Retrieve the pre-calculated per-position divergence values
	reg_missing_masked_sets = get_reg_missing_masked(regions,input_pos_div_file,pop_indexes)
	# stat.write(stat_line + "\n")

	fine_div_win_list=[]
	fine_win_starts=[]
	# Determine whether total number of sites or filtered number of sites will be used to set window size
	if (win_type == 'filtered'):
		# Calculate per-site divergence within windows within regions, gives list of lists of per-window divergence in each region
		fine_div_win_list,fine_win_starts,stat_line = calc_filt_wind_per_region(regions,win_size,reg_missing_masked_sets)
		stat.write(stat_line + "\n")
	if (win_type == 'total'):
		fine_div_win_list,fine_win_starts,stat_line = calc_total_wind_per_region(regions,win_size,reg_missing_masked_sets)
		stat.write(stat_line + "\n")

	print(fine_div_win_list[0:5])
	print(fine_win_starts[0:5])

	return
		

if __name__ == '__main__':
	main()
