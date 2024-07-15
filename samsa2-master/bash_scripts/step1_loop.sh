#! /bin/bash

######################### define flags and usage of the script ####################
function usage() {
    cat <<USAGE

    Usage:
      $0 [-p path to where the bash_scripts directory is found] [-i input directory] [-o output directory]

    Options:
        -p, --path:            Path to where the bash_scripts directory is found
        -i, --input_dir        input directory - path for your scratch folder input_dir
        -o, --output_dir       output_dir - name of the directory you are going to give for the resulting scripts
        -e, --email            your UBC email
        -m, --min_length       minimum_length for TRIMMOMATIC
        -c, --path_common      Path to the common.sh script needed to run the pipeline


USAGE
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

path_bash=
input=
output=
UBC_email=
min_len=
common_script=

while [ "$1" != "" ]; do
    case $1 in
    -p | --path)
        shift
        path_bash=$1
        ;;
		-i | --input_dir)
				shift
				input=$1
				;;
    -o | --output_dir)
    		shift
    		output=$1
    		;;
    -e | --email)
    		shift
    		UBC_email=$1
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

########################## Create script ##############################################

mkdir ${output}

for FILE in $(ls ${input}/*.fastq* | sed 's/_R1*2*.*//' | uniq); do
  rm ${output}/1_$(basename ${FILE}).sh
  touch ${output}/1_$(basename ${FILE}).sh
  cat >> ${output}/1_$(basename ${FILE}).sh << EOL
			#!/bin/bash

			#SBATCH --time=12:00:00
            #SBATCH --nodes=1
            #SBATCH --ntasks=1
            #SBATCH --cpus-per-task=4
            #SBATCH --mem=24G 
            #SBATCH --job-name=step1_$(basename ${FILE})
            #SBATCH --output=o.%x_%A_%a.txt    
            #SBATCH --error=e.%x_%A_%a.txt      
            #SBATCH --mail-user=${UBC_email}
            #SBATCH --mail-type=ALL

			cd ${path_bash};

			export SAMPLE=$(basename $FILE);

			bash step1_final.sh -i ${input}/.. -m ${min_len} -c ${common_script} &> ${input}/../output/step1_$(basename ${FILE}).log
EOL
done


# for FILE in $(ls ${input}/*.fastq* | sed 's/_R1*2*.*//' | uniq); do
# 	echo $(basename ${FILE})
# done
