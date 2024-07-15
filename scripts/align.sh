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

######################### define flags and usage of the script ####################
function usage() {
    cat <<USAGE

    Usage:
      $0 [-p path to analysis directory] [-i input file] [-d input directory] [-db diamond database] [-a aggregation database]

    Options:
        -i, --input        input directory - path for your scratch folder input directory
        -d, --input_directory       Path to where the diamond db is found
        -db, --diamond_db
        -a, --aggregation_db   Path to where your downloaded protein file is found (.faa file)
        -c, --path_common      Path to the common.sh script needed to run the pipeline


USAGE
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

PROJECT=$(pwd)
INPUT_DIR=${PROJECT}/input
input=
diamond_DB=
aggregation_DB=
common_script=

while [ "$1" != "" ]; do
    case $1 in
    -p | --project_directory)
        shift
        PROJECT=$1
        ;;
		-i | --input)
				shift
				input=$1
				;;
    -d | --input_directory)
        shift
        INPUT_DIR=$1
        ;;   
    -db | --diamond_db)
    		shift
    		diamond_DB=$1
    		;;
    -a | --aggregation_db)
    		shift
    		aggregation_DB=$1
    		;;
    -c | --path_common)
    		shift
    		common_script=$1
    		;;
    -h | --help)
        usage
        ;;
    *)
        usage
        exit 1
        ;;
    esac
    shift
done

#############################################

#######################################################################
#
echo -e "NOTE: Before running this script, please run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
#
source "${common_script}"

# actiavte conda env
source ~/.bashrc
conda activate SAMSA2

# Navigate to the project directory
cd $PROJECT 

# Define directory paths
OUT_DIR=${PROJECT}/output
STEP_1="$OUT_DIR/step_1_output"
STEP_2="$OUT_DIR/step_2_output"
STEP_3="$OUT_DIR/step_3_output"
STEP_4="$OUT_DIR/step_4_output"
STEP_5="$OUT_DIR/step_5_output"

# Define number of threads to use. This is automatically set to the number of available CPUs in the Sockeye node
threads=$(nproc)

# Redefine the input file and the file prefix
SAMPLE=${input}
FILENAME=$(basename $SAMPLE | awk -F "_R1" '{print $1}')

# Define database paths
diamond_db=${diamond_DB}
aggregation_db=${aggregation_DB}

####################################################################
# STEP4: ANNOTATING WITH DIAMOND
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

Step=$(grep "ANNOT" ${PROJECT}/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "ANNOT" ]; then

  echo "Now starting on DIAMOND org annotations at: "; date

  for file in $STEP_3/${FILENAME}*ribodepleted.fastq
  do
      shortname=`echo $file | awk -F "ribodepleted" '{print $1 "annotated"}'`
      echo "Now starting on " $file
      echo "Converting to " $shortname

      time checked diamond blastx \
        --db $diamond_db \
        -q $file \
        -a $file.annotated \
        --tmpdir /dev/shm \
        -k 1

      checked diamond view \
      --daa $file.annotated.daa \
      -o $shortname \
      -f tab
  done

  $MKDIR $STEP_4/daa_binary_files

  mv $STEP_3/${FILENAME}*annotated* $STEP_4
  mv $STEP_4/${FILENAME}*.daa $STEP_4/daa_binary_files

  echo "DIAMOND annotations completed at: "; date

  printf "ANNOT\n" >>${PROJECT}/checkpoints/${FILENAME}_checkpoints

else
  printf  "\tThe variable ANNOT is in the checkpoint file. STEP 4 will be skipped.\n"
fi

####################################################################
# STEP5: AGGREGATING WITH ANALYSIS_COUNTER
Step=$(grep "AGGREG" ${PROJECT}/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "AGGREG" ]; then

  for file in $STEP_4/${FILENAME}*annotated
  do
    time checked python $PY_DIR/DIAMOND_analysis_counter.py \
      -I $file \
      -D $aggregation_db \
      -O
    checked python $PY_DIR/DIAMOND_analysis_counter.py \
      -I $file \
      -D $aggregation_db \
      -F
  done

  $MKDIR $STEP_5/Aggreg_results/org_results
  $MKDIR $STEP_5/Aggreg_results/func_results
  mv $STEP_4/${FILENAME}*organism.tsv $STEP_5/Aggreg_results/org_results
  mv $STEP_4/${FILENAME}*function.tsv $STEP_5/Aggreg_results/func_results

  printf "AGGREG\n" >>${PROJECT}/checkpoints/${FILENAME}_checkpoints

else
  printf  "\tThe variable AGGREG is in the checkpoint file. STEP 5 will be skipped.\n"
fi

mkdir -p logs
mv *.logs logs