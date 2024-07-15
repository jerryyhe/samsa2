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
#   1.   FastQC
#
# NOTE: BEFORE running this script, please run package_installation.bash
# and full_database_download.bash located at:
# https://github.com/transcript/samsa2/tree/master/setup in order to set
# up SAMSA2 dependencies and download full databases.
#
#######################################################################
#

######################### define flags and usage of the script ####################
function usage() {
    cat <<USAGE

    Usage:
      $0 [-p path to where the bash_scripts directory is found] [-i input directory] [-o output directory]

    Options:
        -i, --input_dir        input directory - path for your scratch folder input directory
        -m, --min_length       Minimum read length. If below this length, it will be dropped by TRIMMOMATIC
        -c, --path_common      Path to the common.sh script needed to run the pipeline


USAGE
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

input=
min_len=
common_script=

while [ "$1" != "" ]; do
    case $1 in
		-i | --input_dir)
				shift
				input=$1
				;;
    -m | --min_length)
    		shift
    		min_len=$1
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

echo -e "NOTE: Before running this script, please run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
#
source "${common_script}"


# actiavte conda env
module load miniconda3/4.6.14
source activate metatranscriptomics_env

# CURRENT_PROJECT=${input}
CURRENT_PROJECT=${input}

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


####################################################################
#STEP 0.1: create/read checkpoint

mkdir $CURRENT_PROJECT/checkpoints

printf "\nStep 0.1: Checking for the presence of the checkpoint file.\n"
if [ ! -f "$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints" ]
  then
    printf "\tThe file '${CURRENT_SAMPLE}_checkpoints' does not exist in the input directory, creating...\n"
    touch "$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints"
else
    printf "\tThe file '${CURRENT_SAMPLE}_checkpoints' already exists in the input directory.\n"
fi

####################################################################
#STEP 0.2: CHECKING FILES WITH FastQC

mkdir $OUT_DIR/Pre_FastQC

Step=$(grep "PRE_FASTQC" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "PRE_FASTQC" ]
  then

	if ls $INPUT_DIR/${CURRENT_SAMPLE}* &>/dev/null; then
		  for file in $INPUT_DIR/${CURRENT_SAMPLE}*
		    do
                fastqc $file -o $OUT_DIR/Pre_FastQC -t $threads
        done
	fi
printf "PRE_FASTQC\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints
else
  printf  "\tThe variable FASTQC is in the checkpoint file. STEP 1 will be skipped.\n"
fi

####################################################################

# step 0.3 unzipping files with gunzip
#   Note the code for unzipping was not working propperly
#   so we are doing it 'manually'

if ls $INPUT_DIR/${CURRENT_SAMPLE}* &>/dev/null; then
	  for file in $INPUT_DIR/${CURRENT_SAMPLE}*
	    do
              gunzip $file
      done
fi

####################################################################
# STEP 1: CLEANING FILES WITH TRIMMOMATIC

mkdir $OUT_DIR/step_1_output

Step=$(grep "TRIMMO" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "TRIMMO" ]
  then

		if ls $INPUT_DIR/${CURRENT_SAMPLE}*.fastq &>/dev/null; then
		  for file in $INPUT_DIR/${CURRENT_SAMPLE}*.fastq
		  do
			gunzip $file
		  done
		fi

		$MKDIR $STEP_1
		paired=false
		for f in $INPUT_DIR/${CURRENT_SAMPLE}*R1*q
		do
			f2=`echo $f | awk -F "R1" '{print $1 "R2" $2}'`
    		out_path=`echo $f | awk -F "_R1" '{print $1 ".cleaned"}'`
			if [ -f $f2 ]; then
			  paired=true
			  time java -jar $TRIMMOMATIC PE -phred33 -threads $threads $f $f2 \
				$out_path".forward" $out_path".forward_unpaired" $out_path".reverse" $out_path".reverse_unpaired" \
				ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:${min_len}
			else
			  time java -jar $TRIMMOMATIC SE -phred33 -threads $threads $f $out_path SLIDINGWINDOW:4:15 MINLEN:${min_len}
			fi
		done

		if $paired; then
		  mv "$INPUT_DIR/${CURRENT_SAMPLE}.cleaned.forward"* $OUT_DIR/step_1_output
		  mv "$INPUT_DIR/${CURRENT_SAMPLE}.cleaned.reverse"* $OUT_DIR/step_1_output
		else
		  mv "$INPUT_DIR/${CURRENT_SAMPLE}.cleaned" $OUT_DIR/step_1_output
		fi
		printf "TRIMMO\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints

else
  printf  "\tThe variable TRIMMO is in the checkpoint file. STEP 1 will be skipped.\n"
fi

####################################################################

#STEP 1.1: CHECKING FILES WITH FastQC post Trimmomatic

mkdir $OUT_DIR/Post_FastQC

Step=$(grep "POST_FASTQC" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "POST_FASTQC" ]
  then
	if ls $OUT_DIR/step_1_output/${CURRENT_SAMPLE}*.cleaned* &>/dev/null; then
		  for file in $OUT_DIR/step_1_output/${CURRENT_SAMPLE}*.cleaned*
		    do
                fastqc $file -o $OUT_DIR/Post_FastQC -t $threads
        done
	fi
printf "POST_FASTQC\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints
else
  printf  "\tThe variable Post_FASTQC is in the checkpoint file. STEP 1 will be skipped.\n"
fi

####################################################################

#
# STEP 2: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#       Example: control_1.R1.fastq
#                control_1.R2.fastq

Step=$(grep "MERGING" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "MERGING" ]
  then

$MKDIR $STEP_2
if $paired; then
  for file in $STEP_1/${CURRENT_SAMPLE}*.cleaned.forward
  do
    f2=`echo $file | awk -F "cleaned.forward" '{print $1 "cleaned.reverse"}'`
    shortname=`echo $file | awk -F "cleaned.forward" '{print $1 "merged"}'`
    time $PEAR -f $file -r $f2 -j $threads -o $STEP_2/${shortname##*/}
  done
else
  for file in $STEP_1/${CURRENT_SAMPLE}*.cleaned
  do
    new_name=`echo $file | awk -F "cleaned" '{print $1 "merged.assembled.fastq"}'`
    cp $file $STEP_2/${new_name##*/}
  done
fi

printf "MERGING\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints
else
  printf  "\tThe variable MERGING is in the checkpoint file. STEP 2 will be skipped.\n"
fi

####################################################################

# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.
Step=$(grep "RAW" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "RAW" ]
  then

if [[ -f $STEP_2/raw_counts.txt ]]; then
    rm $STEP_2/raw_counts.txt
fi
touch $STEP_2/raw_counts.txt

if $paired; then
  for file in $STEP_1/${CURRENT_SAMPLE}*cleaned.forward
  do
    time python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
  done
else
  for file in $STEP_1/${CURRENT_SAMPLE}*cleaned
  do
    time python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
  done
fi

printf "RAW\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints

else
  printf  "\tThe variable RAW is in the checkpoint file. STEP 2.9 will be skipped.\n"
fi

####################################################################

# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).
Step=$(grep "RIBO" $CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints)
if [ "${Step}" != "RIBO" ]
  then

for file in $STEP_2/${CURRENT_SAMPLE}*.assembled.fastq
do
  shortname=`echo $file | awk -F "assembled" '{print $1 "ribodepleted"}'`
  time $SORTMERNA -a $threads \
    --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db \
    --reads $file --aligned $file.ribosomes --other $shortname --fastx \
    --log -v -m 4096
done

$MKDIR $STEP_3
mv $STEP_2/${CURRENT_SAMPLE}*ribodepleted* $STEP_3

printf "RIBO\n" >>$CURRENT_PROJECT/checkpoints/${CURRENT_SAMPLE}_checkpoints

else
  printf  "\tThe variable RIBO is in the checkpoint file. STEP 3 will be skipped.\n"
fi
