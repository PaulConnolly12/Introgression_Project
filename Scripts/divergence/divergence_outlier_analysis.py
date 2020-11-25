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
		for (lnum, line) in enumerate(mel_file):
			div_vals = line.strip().split('\t')
			for i in pops_of_interest:




	return


# Determine regions to search through at a finer scale
def find_outliers(in_window_div,in_window_starts):

	return regions


# Calculate filtered windows across each specified region
def calc_filt_wind_per_region(regions,win_size):

	return



# Input arguments
def parse_args():

	# Options for input files
	parser = argparse.ArgumentParser()

	parser.add_argument('--window_pop_list', nargs='+', required=False,
						help='Provide the list of populations used in generating the ')
	parser.add_argument('--window_divergence_file', nargs='+', required=True,
						help='Provide the window divergence file to be used in analysis')
	parser.add_argument('--window_starts_file', nargs='+', required=True,
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
	input_window_file = str(args.window_starts_file)
	input_start_file = str(args.window_starts_file)
	out = str(args.out) + "_" + str(args.window_type) + "_" + str(args.window_size) + ".fine_div_window"
	start = str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".win_starts"
	stat = open((str(args.out)  + "_" + str(args.window_type) + "_" + str(args.window_size) + ".status"), 'w')
	win_type = str(args.window_type)
	size = int(args.window_size)

	stat.write("Running a fine-scale divergence analysis at " + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + "\n")
	stat.write("Window size: " + str(size) + "\n" + "Window type: " + str(win_type) + "\n")

	# Manually define the indexes of the population pair of interest
	pops_of_interest = [1,21]

	# Parse the inputs
	in_window_div = parse_win_div(input_window_file,pops_of_interest)
	in_window_starts = parse_win_starts(input_start_file)

	# Manually define the proportion of the windows to pick as outliers
	outlier_proportion = 0.025

	# Identify all outlier regions
	pop_list,stat_line = find_outliers(in_window_div,in_window_starts)
	stat.write(stat_line + "\n")

	# Determine whether total number of sites or filtered number of sites will be used to set window size
	if (win_type == 'total'):
		# Calculate per-site divergence within windows within regions, gives list of lists of per-window divergence in each region
		
	if (win_type == 'filtered'):
		

if __name__ == '__main__':
	main()
