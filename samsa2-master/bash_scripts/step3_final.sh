#!/bin/bash

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# master_script.sh
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified March 12, 2019
#
####################################################################
#
# This script sets up and runs through ALL steps in the SAMSA pipeline
# before the analysis (which is done in R, likely in RStudio).  Each
# step is set up below.
#
# The steps are:
#   1.   Read cleaning with Trimmomatic
#   2.   Merging with PEAR, if applicable
#   3.   rRNA removal with SortMeRNA
#   4.   Annotation using DIAMOND (by default against the RefSeq database)
#   5.   Aggregation using analysis_counter.py
#   4.1  Annotation using DIAMOND against the Subsystems database
#   5.1  Aggregation using Subsystems-specific analysis counter.py
#   6.   Running R scripts to get DESeq statistical analysis.
#
# NOTE: BEFORE running this script, please run package_installation.bash
# and full_database_download.bash located at:
# https://github.com/transcript/samsa2/tree/master/setup in order to set
# up SAMSA2 dependencies and download full databases.
#
#######################################################################
#
echo -e "NOTE: Before running this script, please run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
#
source "/home/dmpepin/samsa2-master/bash_scripts/lib/common.sh"


CURRENT_PROJECT=/scratch/st-ctropini-1/deannabigtest

INPUT_DIR=${CURRENT_PROJECT}/input
OUT_DIR=${CURRENT_PROJECT}/output

# number of threads
threads=$NCPUS

STEP_1="$OUT_DIR/step_1_output"
STEP_2="$OUT_DIR/step_2_output"
STEP_3="$OUT_DIR/step_3_output"
STEP_4="$OUT_DIR/step_4_output"
STEP_5="$OUT_DIR/step_5_output"

CURRENT_SAMPLE=$SAMPLE

##################################################################
#
# At this point, all the results files are ready for analysis using R.
# This next step performs basic DESeq2 analysis of the RefSeq organism, function,
# and Subsystems annotations.
#
# More complex R analyses may be performed using specific .sh analysis scripts.
#
# STEP 6: R ANALYSIS
# Note: For R to properly identify files to compare/contrast, they must include
# the appropriate prefix (either "control_$file" or experimental_$file")!
Step=$(grep "R_ANALYSIS" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "R_ANALYSIS" ]
  then

checked Rscript $R_DIR/run_DESeq_stats.R \
  -I $STEP_5/RefSeq_results/org_results \
  -O RefSeq_org_DESeq_results.tab \
  -R $STEP_2/raw_counts.txt
checked Rscript $R_DIR/run_DESeq_stats.R \
  -I $STEP_5/RefSeq_results/func_results \
  -O RefSeq_func_DESeq_results.tab \
  -R $STEP_2/raw_counts.txt
checked Rscript $R_DIR/Subsystems_DESeq_stats.R \
  -I $STEP_5/Subsystems_results \
  -O Subsystems_level-1_DESeq_results.tab -L 1 \
  -R $STEP_2/raw_counts.txt

printf "R_ANALYSIS\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints

else
  printf  "\tThe variable R_ANALYSIS is in the checkpoint file. STEP 6 will be skipped.\n"
fi

echo "Master bash script finished running at: "; date
exit
