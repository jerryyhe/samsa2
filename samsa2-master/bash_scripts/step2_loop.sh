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
        -d, --diamond_db       Path to where the diamond db is found
        -a, --aggregation_db   Path to where your downloaded protein file is found (.faa file)
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
diamond_DB=
aggregation_DB=
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
    -d | --diamond_db)
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

########################## Create script ##############################################

mkdir ${output}

for FILE in $(ls ${input}/*.fastq* | sed 's/_R1*2*.*//' | uniq); do
  rm ${output}/2_$(basename ${FILE}).sh
  touch ${output}/2_$(basename ${FILE}).sh
  cat >> ${output}/2_$(basename ${FILE}).sh << EOL
			#!/bin/bash

            #SBATCH --time=5:00:00
            #SBATCH --nodes=1
            #SBATCH --ntasks=1
            #SBATCH --cpus-per-task=4
            #SBATCH --mem=64G 
            #SBATCH --job-name=step2_$(basename ${FILE})
            #SBATCH --output=o.%x_%A_%a.txt    
            #SBATCH --error=e.%x_%A_%a.txt      
            #SBATCH --mail-user=${UBC_email}
            #SBATCH --mail-type=ALL

			cd ${path_bash};

			export SAMPLE=$(basename $FILE);

			bash step2_final.sh -i ${input}/.. -d ${diamond_DB} -a ${aggregation_DB} -c ${common_script} &> ${input}/../output/step2_$(basename ${FILE}).log
EOL
done
