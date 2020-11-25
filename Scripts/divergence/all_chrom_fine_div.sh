# Purpose: calculate per-site divergence from D. sim for arbitrary number of D. mel populations,
# across a set of regions of interest using the divergence.py python script

# echo dirname "$0"


# Run an instance with the following command:
# sh all_mel_sim_divergence.sh
## python divergence.py --d_mel_count_files <list of files> --d_sim_ref_count_file <file> --d_mel_count_file_path <path to files> --d_sim_count_file_path <path to file> --window_size <int> --window_type <filtered or total> --out <path + filename>

# Define the parameters, in case the script will be modified later to take in a specific custom parameter
POP_COUNT_PREFIXES=("AllCounts_CO" "AllCounts_ZI")
# COUNT_FILE_DIR="/Users/cmcallester/Documents/Pool Lab/Introgression/Introgression_Project/Data/Counts/"
COUNT_FILE_DIR="../../Data/Counts/"
SIM_REF_PREFIX="AllCounts_sim"
SIM_DIR=$COUNT_FILE_DIR
OUT_FILE_PREFIX="fine_div_test"
OUT_DIR="Outputs/fine_scale_test/"
WIN_SIZE=100
WIN_TYPE="filtered"
# WORKING_DIR="/Users/cmcallester/Documents/Pool Lab/Introgression/Introgression_Project/"
WORKING_DIR="../../"
DIV_ANALYSIS_SCRIPT="Scripts/divergence/divergence_outlier_analysis.py"

# INPUT_DIR="/Users/cmcallester/Documents/Pool Lab/Introgression/Introgression_Project/Data/divergence/"
INPUT_DIR="../../Data/divergence/"
INPUT_PREFIX="all_pop_new"
INPUT_WIN_TYPE="filtered"
INPUT_WIN_SIZE=100000


# Take a time-stamp for log recording
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# Run each chrom:
for CHROM in X 2L 2R 3L 3R
do
	# Prepare chromosome-specific pop file names
	POP_COUNT_FILES=()
	for i in ${!POP_COUNT_PREFIXES[@]}; do
		POP_FILE=${POP_COUNT_PREFIXES[i]}'_Chr'$CHROM'.txt'
		POP_COUNT_FILES+=($POP_FILE)
	done
	SIM_REF=$SIM_REF_PREFIX'_Chr'$CHROM'.txt'
	# Prepare output file names
	OUT=$WORKING_DIR$OUT_DIR$OUT_FILE_PREFIX'_Chr'$CHROM
	NOHUP_FILE='nohup_Chr'$CHROM'_'$TIMESTAMP'.out'
	# Prepare input file names
	INPUT_DIV_FILE=$INPUT_DIR$INPUT_PREFIX'_Chr'$CHROM'_'$INPUT_WIN_TYPE'_'$INPUT_WIN_SIZE'.div_window'
	INPUT_STARTS_FILE=$INPUT_DIR$INPUT_PREFIX'_Chr'$CHROM'_'$INPUT_WIN_TYPE'_'$INPUT_WIN_SIZE'.win_starts'

	# echo python "'"$WORKING_DIR$DIV_ANALYSIS_SCRIPT"'" --window_divergence_file "'"$INPUT_DIV_FILE"'" --window_starts_file "'"$INPUT_STARTS_FILE"'" --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path "'"$COUNT_FILE_DIR"'" --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT
	echo python $WORKING_DIR$DIV_ANALYSIS_SCRIPT --window_divergence_file $INPUT_DIV_FILE --window_starts_file $INPUT_STARTS_FILE --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path "'"$COUNT_FILE_DIR"'" --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT
	python $WORKING_DIR$DIV_ANALYSIS_SCRIPT --window_divergence_file $INPUT_DIV_FILE --window_starts_file $INPUT_STARTS_FILE --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT
	# echo nohup python $WORKING_DIR$DIV_ANALYSIS_SCRIPT --window_divergence_file $INPUT_DIV_FILE --window_starts_file $INPUT_STARTS_FILE --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT ">" $NOHUP_FILE "2>&1"
	# nohup python $WORKING_DIR$DIV_ANALYSIS_SCRIPT --window_divergence_file $INPUT_DIV_FILE --window_starts_file $INPUT_STARTS_FILE --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT > $NOHUP_FILE 2>&1 &
done
