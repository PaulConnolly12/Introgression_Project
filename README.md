Simulans Introgression Project Information
---------------------------------------------------------------------------------------------

## A repository for scripting associated with analyzing differential divergence of cosmopolital Drosophila melanogaster populations to the ancestral Zambian population. UNDER ACTIVE DEVELOPMENT

This document describes the location and usage of relevant scripts as well as their output

## Directory Information on GenePool Lab Server

Working Directory: /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/
Scripts: /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/scripts/
Test Data: /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/test_data/
Full Data: /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/data/
Output and Results : /Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/EKH/results/

"Call me Ishmael. Some years ago—never mind how long precisely—having little or no money in my purse, and nothing particular to interest me on shore, I thought I would sail about a little and see the watery part of the world. It is a way I have of driving off the spleen and regulating the circulation. Whenever I find myself growing grim about the mouth; whenever it is a damp, drizzly November in my soul; whenever I find myself involuntarily pausing before coffin warehouses, and bringing up the rear of every funeral I meet; and especially whenever my hypos get such an upper hand of me, that it requires a strong moral principle to prevent me from deliberately stepping into the street, and methodically knocking people’s hats off—then, I account it high time to get to sea as soon as I can. This is my substitute for pistol and ball. With a philosophical flourish Cato throws himself upon his sword; I quietly take to the ship. There is nothing surprising in this. If they but knew it, almost all men in their degree, some time or other, cherish very nearly the same feelings towards the ocean with me."

Contributors
----------------------------------
Emma Howell\
Paul Connoly\
Chris McAllester\


Available Scripts
----------------------------------

There are four scripts relevant to this project, two python scripts to calculate divergence and perform the ABBA BABA test as well as two R scripts for plotting the divergence results and testing for a correlation between divergence and recombination

### Divergence Script

Filename: divergence.py

Summary: This script takes in a variable number of D. melanogaster count files as well as a D. simulans reference count file and calculates per-site divergence across the full chromosome arm as well as on a per-window basis. Outputs window start locations. Count files for input should be organized as a line per site of tab-separated counts of bases as 'A T G C N'

Usage: Details about running script can be found in script header. Script can be run from python, using the test data, as:

	python Scripts\divergence\divergence.py --d_mel_count_files testdata_CO_Chr2R.txt testdata_FR_Chr2R.txt --d_sim_ref_count_file testdata_sim_Chr2R.txt --d_mel_count_file_path test_data\ --d_sim_count_file_path test_data\ --window_size 10 --window_type total --out test_outputs\test_run 


### ABBA BABA D-statistic Script

Filename: abba_baba.py

Summary: This script takes in count files for a P1, P2, P3, and PO population then calculates Patterson's D statistic according to both a count-based and frequency-based approach

Usage: Details about running script can be found in script header

Notes: In order to evaluate the significance of any D value greater than 0, a blcok jackknife procedure can be performed. This would require calculating D in particular windows and seeing if the calculated D value is robust to changing which windows are included. Currently, the script is not equipped to calulcate D on a window basis but there are comments within the script about how a window-based approach could be implemented

### R Scripts

Filename: sims_div_plots.R and introgression_correlation.R

Summary: These scripts contain the code used to make plots from the divergence output and to perform the correlation test using the recombination rate data

Usage: These scripts are not currently set up to be run on the command line and have filenames and paths hard-coded. Instead, they need to be downloaded and manipulated using R or Rstudio. This would also require transferring the input files from genepool onto whichever computer is being used to modify the scripts

Notes: These scripts provide the code used to make plots and perform the correlation test but are not meant to be run in their entirety. Instead, blocks of code from these scripts can be re-purposed for general plotting purposes


Data Information
----------------------------------

### Test Data

These representative count files each contain 100 lines and are constructed from their full counterparts. These files are useful for testing python scripts that involve performing calculations on the full count files. The small size of these files allow you to print relevant output to the screen to ensure your code is working as expected

### Data

While some of the count files have been copied to this directory redundantly, this folder provides easy access to the processed yakuba count files and contains the recombination rate data from Comeron (2012) 

Results Information
----------------------------------

### Divergence

This sub-directory contains the output from the divergence script that has been run for all populations; for chromosomes Chr2R, Chr2L, Chr3R, Chr3L, ChrX; for window sizes of 100kb and 1kb; and for window options of filtered and total

### Patterson's D

This sub-directory contains the output from the ABBA BABA script that has been run for chromosomes Chr2R, Chr2L, Chr3R, Chr3L, ChrX under a scenario where P1 is CO, P2 is ZI, P3 is D. simulans, and PO is D. yakuba

### Plots

This sub-directory contains all of the plots created using the above output data and corresponding R scripts 

Miscellaneous Information
----------------------------------

The following commands were used to create the toy data files

```bash
cat AllCounts_FR_Chr2R.txt | awk '{if(($1 != 0) || ($2 != 0) || ($3 != 0) || ($4 != 0)) print $1"\t"$2"\t"$3"\t"$4}' | head -n100 > ../EKH/test_data/toydata_FR_Chr2R.txt
cat AllCounts_CO_Chr2R.txt | awk '{if(($1 != 0) || ($2 != 0) || ($3 != 0) || ($4 != 0)) print $1"\t"$2"\t"$3"\t"$4}' | head -n100 > ../EKH/test_data/toydata_CO_Chr2R.txt
cat AllCounts_sim_Chr2R.txt | awk '{if(($1 != 0) || ($2 != 0) || ($3 != 0) || ($4 != 0)) print $1"\t"$2"\t"$3"\t"$4}' | head -n100 > ../EKH/test_data/toydata_sim_Chr2R.txt
cat AllCounts_yakuba_Chr2R.txt | awk '{if(($1 != 0) || ($2 != 0) || ($3 != 0) || ($4 != 0)) print $1"\t"$2"\t"$3"\t"$4}' | head -n100 > ../EKH/test_data/toydata_yak_Chr2R.txt
```


Installation
----------------------------------

Download the scripts from this repository:

	git clone https://github.com/MesserLab/SLiM.git

Scripts present require python or R to run. 
