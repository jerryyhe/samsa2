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

######################### Define Flags and Script Help #########################
function usage() {
    cat <<USAGE

    Usage:
      $0 [-p path to analysis directory] [-d input directory] [-i input file]

    Options:
        -p, --project_directory   Path of the directory in which you will perform this analysis. Default is the current working directory.
        -d, --input_directory     Input directory. Path to the location of your input reads. Ensure that this directory ONLY contains your input reads and no other files. Default is a directory called input in your project directory.
        -i, --input               Input FASTQ file. 
        -m, --min_length          The minimum read length used for TRIMMOMATIC. If below this length, it will be dropped by TRIMMOMATIC
        -c, --path_common         Path to the common.sh script needed to run the pipeline. Default is _____


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
min_len=
common_script=

while [ "$1" != "" ]; do
    case $1 in
    -p | --project_directory)
				shift
				PROJECT=$1
				;;
    -d | --input_directory)
        shift
        INPUT_DIR=$1
        ;;   
		-i | --input)
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

# VARIABLES - set starting location and starting files location pathways
source "${common_script}"

# Activate Conda Environment
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

#- Note that section headers with a '#-' were modified from the original SAMSA2 code
####################################################################
#STEP 0.1: CREATE CHECKPOINTS

mkdir -p $PROJECT/checkpoints

printf "\nStep 0.1: Checking for the presence of the checkpoint file.\n"
if [ ! -f "$PROJECT/checkpoints/${FILENAME}_checkpoints" ]; then
    printf "\tThe file '${FILENAME}_checkpoints' does not exist in the input directory, creating...\n"
    touch "$PROJECT/checkpoints/${FILENAME}_checkpoints"
else
    printf "\tThe file '${FILENAME}_checkpoints' already exists in the input directory.\n"
fi

####################################################################
#STEP 0.2: CHECKING FILES WITH FastQC

mkdir -p $OUT_DIR/Pre_FastQC

Step=$(grep "PRE_FASTQC" $PROJECT/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "PRE_FASTQC" ]; then

	if ls $INPUT_DIR/${FILENAME}* &>/dev/null; then
		  for file in $INPUT_DIR/${FILENAME}*
		    do
                fastqc $file -o $OUT_DIR/Pre_FastQC -t $threads
        done
	fi
printf "PRE_FASTQC\n" >>$PROJECT/checkpoints/${FILENAME}_checkpoints
else
  printf  "\tThe variable FASTQC is in the checkpoint file. STEP 1 will be skipped.\n"
fi

####################################################################
#-STEP 0.3 UNZIPPING GZIPPED FILES
#   Note the code for unzipping was not working properly so we are doing it 'manually'

if file ${SAMPLE} | grep -q "gzip compressed data"; then
  gunzip $SAMPLE
  gunzip $(echo $SAMPLE | awk -F "R1" '{print $1 "R2" $2}')
fi

####################################################################
#-STEP 1: CLEANING FILES WITH TRIMMOMATIC

mkdir -p $OUT_DIR/step_1_output

Step=$(grep "TRIMMO" $PROJECT/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "TRIMMO" ]; then

		$MKDIR $STEP_1
		paired=false
		for f in ${INPUT_DIR}/${FILENAME}*; do

			f2=$(echo $f | awk -F "R1" '{print $1 "R2" $2}')
      out_path=$(echo $f | awk -F "_R1" '{print $1 ".cleaned"}')

			if [ -f $f2 ]; then
			  paired=true
			  time java -jar $TRIMMOMATIC PE \
        -phred33 \
        -threads $threads \
        $f \
        $f2 \
				$out_path".forward" \
        $out_path".forward_unpaired" \
        $out_path".reverse" \
        $out_path".reverse_unpaired" \
				ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:15 \
        MINLEN:${min_len}

			else
			  time java -jar $TRIMMOMATIC SE \
        -phred33 \
        -threads $threads \
        $f \
        $out_path \
        SLIDINGWINDOW:4:15 \
        MINLEN:${min_len}

			fi
		done

		if $paired; then
		  mv "$INPUT_DIR/${FILENAME}.cleaned.forward"* $OUT_DIR/step_1_output
		  mv "$INPUT_DIR/${FILENAME}.cleaned.reverse"* $OUT_DIR/step_1_output
		else
		  mv "$INPUT_DIR/${FILENAME}.cleaned" $OUT_DIR/step_1_output
		fi
		printf "TRIMMO\n" >>$PROJECT/checkpoints/${FILENAME}_checkpoints

else
  printf  "\tThe variable TRIMMO is in the checkpoint file. STEP 1 will be skipped.\n"
fi

####################################################################
#-STEP 1.1: CHECKING FILES WITH FastQC post Trimmomatic

mkdir -p $OUT_DIR/Post_FastQC

Step=$(grep "POST_FASTQC" $PROJECT/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "POST_FASTQC" ]
  then
	if ls $OUT_DIR/step_1_output/${FILENAME}*.cleaned* &>/dev/null; then
		  for file in $OUT_DIR/step_1_output/${FILENAME}*.cleaned*
		    do
                fastqc $file -o $OUT_DIR/Post_FastQC -t $threads
        done
	fi
printf "POST_FASTQC\n" >>$PROJECT/checkpoints/${FILENAME}_checkpoints
else
  printf  "\tThe variable Post_FASTQC is in the checkpoint file. STEP 1 will be skipped.\n"
fi

####################################################################
# STEP2: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#       Example: control_1_R1.fastq
#                control_1_R2.fastq

Step=$(grep "MERGING" $PROJECT/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "MERGING" ]; then

  $MKDIR $STEP_2
  if $paired; then
    for file in $STEP_1/${FILENAME}*.cleaned.forward
    do
      f2=`echo $file | awk -F "cleaned.forward" '{print $1 "cleaned.reverse"}'`
      shortname=`echo $file | awk -F "cleaned.forward" '{print $1 "merged"}'`
      time $PEAR -f $file -r $f2 -j $threads -o $STEP_2/${shortname##*/}
    done
  else
    for file in $STEP_1/${FILENAME}*.cleaned
    do
      new_name=`echo $file | awk -F "cleaned" '{print $1 "merged.assembled.fastq"}'`
      cp $file $STEP_2/${new_name##*/}
    done
  fi

  printf "MERGING\n" >>$PROJECT/checkpoints/${FILENAME}_checkpoints

else
  printf  "\tThe variable MERGING is in the checkpoint file. STEP 2 will be skipped.\n"
fi

####################################################################
#STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.
Step=$(grep "RAW" $PROJECT/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "RAW" ]
  then

if [[ -f $STEP_2/raw_counts.txt ]]; then
    rm $STEP_2/raw_counts.txt
fi
touch $STEP_2/raw_counts.txt

if $paired; then
  for file in $STEP_1/${FILENAME}*cleaned.forward
  do
    time python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
  done
else
  for file in $STEP_1/${FILENAME}*cleaned
  do
    time python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
  done
fi

printf "RAW\n" >>$PROJECT/checkpoints/${FILENAME}_checkpoints

else
  printf  "\tThe variable RAW is in the checkpoint file. STEP 2.9 will be skipped.\n"
fi

####################################################################
# STEP3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).
Step=$(grep "RIBO" $PROJECT/checkpoints/${FILENAME}_checkpoints)
if [ "${Step}" != "RIBO" ]
  then

for file in $STEP_2/${FILENAME}*.assembled.fastq
do
  shortname=`echo $file | awk -F "assembled" '{print $1 "ribodepleted"}'`
  time $SORTMERNA -a $threads \
    --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db \
    --reads $file --aligned $file.ribosomes --other $shortname --fastx \
    --log -v -m 4096
done

$MKDIR $STEP_3
mv $STEP_2/${FILENAME}*ribodepleted* $STEP_3

printf "RIBO\n" >>$PROJECT/checkpoints/${FILENAME}_checkpoints

else
  printf  "\tThe variable RIBO is in the checkpoint file. STEP 3 will be skipped.\n"
fi
