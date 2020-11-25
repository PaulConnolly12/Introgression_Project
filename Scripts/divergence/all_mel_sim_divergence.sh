# Purpose: calculate per-site divergence between D. sim for arbitrary number of D. mel populations,
# across each chromosome using the divergence.py python script

# echo dirname "$0"


# Run an instance with the following command:
# sh all_mel_sim_divergence.sh
## python divergence.py --d_mel_count_files <list of files> --d_sim_ref_count_file <file> --d_mel_count_file_path <path to files> --d_sim_count_file_path <path to file> --window_size <int> --window_type <filtered or total> --out <path + filename>

# Define the parameters, in case the script will be modified later to take in a specific custom parameter
POP_COUNT_PREFIXES=("AllCounts_B" "AllCounts_CO" "AllCounts_EA" "AllCounts_EF" "AllCounts_EG" "AllCounts_FR" "AllCounts_GA" "AllCounts_GU" "AllCounts_I" "AllCounts_KF" "AllCounts_Kenya" "AllCounts_MW" "AllCounts_NG" "AllCounts_N" "AllCounts_RAL" "AllCounts_RG" "AllCounts_SD" "AllCounts_SP" "AllCounts_T" "AllCounts_WAF" "AllCounts_W" "AllCounts_ZI")
COUNT_FILE_DIR="/Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/"
SIM_REF_PREFIX="AllCounts_sim"
SIM_DIR=$COUNT_FILE_DIR
OUT_FILE_PREFIX="all_pop_new"
OUT_DIR="New_Outputs/divergence/"
WIN_SIZE=100000
WIN_TYPE="filtered"
WORKING_DIR="/Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/PFC/introgression_project_old/"
DIV_SCRIPT="Scripts/divergence.py"

# Take a time-stamp for log recording
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# Run each chrom:
for CHROM in X 2L 2R 3L 3R
do
	POP_COUNT_FILES=()
	for i in ${!POP_COUNT_PREFIXES[@]}; do
		POP_FILE=${POP_COUNT_PREFIXES[i]}'_Chr'$CHROM'.txt'
		POP_COUNT_FILES+=($POP_FILE)
	done
	OUT=$WORKING_DIR$OUT_DIR$OUT_FILE_PREFIX'_Chr'$CHROM
	SIM_REF=$SIM_REF_PREFIX'_Chr'$CHROM'.txt'
	NOHUP_FILE='nohup_Chr'$CHROM'_'$TIMESTAMP'.out'
	# echo python $WORKING_DIR$DIV_SCRIPT --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT
	# python $WORKING_DIR$DIV_SCRIPT --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT
	echo nohup python $WORKING_DIR$DIV_SCRIPT --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT ">" $NOHUP_FILE "2>&1"
	nohup python $WORKING_DIR$DIV_SCRIPT --d_mel_count_files ${POP_COUNT_FILES[@]} --d_sim_ref_count_file $SIM_REF --d_mel_count_file_path $COUNT_FILE_DIR --d_sim_count_file_path $SIM_DIR --window_size $WIN_SIZE --window_type $WIN_TYPE --out $OUT > $NOHUP_FILE 2>&1 &
done

## python divergence.py --d_mel_count_files AllCounts_B_Chr2R.txt AllCounts_CO_Chr2R.txt AllCounts_EA_Chr2R.txt AllCounts_EF_Chr2R.txt AllCounts_EG_Chr2R.txt AllCounts_FR_Chr2R.txt AllCounts_GA_Chr2R.txt AllCounts_GU_Chr2R.txt AllCounts_I_Chr2R.txt AllCounts_KF_Chr2R.txt AllCounts_Kenya_Chr2R.txt AllCounts_MW_Chr2R.txt AllCounts_NG_Chr2R.txt AllCounts_N_Chr2R.txt AllCounts_RAL_Chr2R.txt AllCounts_RG_Chr2R.txt AllCounts_SD_Chr2R.txt AllCounts_SP_Chr2R.txt AllCounts_T_Chr2R.txt AllCounts_WAF_Chr2R.txt AllCounts_W_Chr2R.txt AllCounts_ZI_Chr2R.txt --d_sim_ref_count_file AllCounts_sim_Chr2R.txt --d_mel_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/ --d_sim_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/counts/ --window_size 10000 --window_type filtered --out test_filtered
## python divergence.py --d_mel_count_files toydata_CO_Chr2R.txt toydata_FR_Chr2R.txt --d_sim_ref_count_file toydata_sim_Chr2R.txt --d_sim_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/ --d_mel_count_file_path /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/ --window_size 10 --window_type total --out /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_run 
