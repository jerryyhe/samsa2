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

# Diamond databases
diamond_database="/project/st-ctropini-1/dpepin/metatranscriptomics_BtS247/full_databases/S24-7_Bt_sep_annotations.dmnd"

# Aggregation databases

aggregation_db="/project/st-ctropini-1/dpepin/metatranscriptomics_BtS247/full_databases/S24-7_Bt_sep_annotations.faa"

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.
Step=$(grep "ANNOT" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "ANNOT" ]
  then

echo "Now starting on DIAMOND org annotations at: "; date

for file in $STEP_3/${CURRENT_SAMPLE}*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "annotated"}'`
    echo "Now starting on " $file
    echo "Converting to " $shortname
    time checked $DIAMOND blastx --db $diamond_database -q $file -a $file.annotated --tmpdir /dev/shm -k 1
    checked $DIAMOND view --daa $file.annotated.daa -o $shortname -f tab
done

$MKDIR $STEP_4/daa_binary_files

mv $STEP_3/${CURRENT_SAMPLE}*annotated* $STEP_4
mv $STEP_4/${CURRENT_SAMPLE}*.daa $STEP_4/daa_binary_files

echo "DIAMOND annotations completed at: "; date

printf "ANNOT\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints

else
  printf  "\tThe variable ANNOT is in the checkpoint file. STEP 4 will be skipped.\n"
fi

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER
Step=$(grep "AGGREG" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "AGGREG" ]
  then

for file in $STEP_4/${CURRENT_SAMPLE}*annotated
do
  time checked python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $aggregation_db -O
  checked python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $aggregation_db -F
done

$MKDIR $STEP_5/Aggreg_results/org_results
$MKDIR $STEP_5/Aggreg_results/func_results
mv $STEP_4/${CURRENT_SAMPLE}*organism.tsv $STEP_5/Aggreg_results/org_results
mv $STEP_4/${CURRENT_SAMPLE}*function.tsv $STEP_5/Aggreg_results/func_results

printf "AGGREG\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints

else
  printf  "\tThe variable AGGREG is in the checkpoint file. STEP 5 will be skipped.\n"
fi
