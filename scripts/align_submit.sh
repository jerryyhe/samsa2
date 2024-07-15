#! /bin/bash

# Define and export the directory of containing the SAMSA2 scripts and databases
export SAMSA2_DIR=$(dirname "$(realpath "$0")")/..

######################### Define Flags and Script Help #########################
function usage() {
    cat <<USAGE

    Usage:
      $0 [-p path to analysis directory] [-i reads input directory] 

    Options:
        -p, --project_directory:            Path of the directory in which you will perform this analysis. Default is the current working directory.
        -i, --input_directory               Input directory. Path to the location of your input reads. Ensure that this directory ONLY contains your input reads and no other files.
        -e, --email                         Your UBC email.
        -d, --diamond_db                    DIAMOND database to use for alignment.
        -a, --aggregation_db                Database to use for read aggregation. This is usually the .faa file used to make a DIAMOND database.
        -c, --path_common                   Path to the common.sh script needed to run the pipeline. Default is _____
        -t, --runtime                       The amount of time requested for read preprocessing in hh:mm:ss format. Default is 5:00:00
        -A, --allocation                    The allocation code of your sockeye allocation owner. Default is st-ctropini-1
        -ncpus                              The number of CPUs desired for read preprocessing. Default is 8
        -mem, --memory                      The amount of RAM desired for read preprocessing. Default is 64.


USAGE
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# Define default variables
PROJECT=$(pwd)
INPUT=
TIME=5:00:00
ALLOCATION=st-ctropini-1
ncpus=8
nram=64
EMAIL=
diamond_DB=
aggregation_DB=
common_script=${SAMSA2_DIR}/samsa2-master/bash_scripts/lib/common.sh

while [ "$1" != "" ]; do
    case $1 in
    -p | --project_directory)
        shift
        PROJECT=$1
        ;;
    -i | --input_directory)
            shift
            INPUT=$1
            ;;
    -e | --email)
    		shift
    		EMAIL=$1
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
    -t | --runtime)
    		shift
    		TIME=$1
    		;;
    -a | --allocation)
            shift
            ALLOCATION=$1
            ;;
    -ncpus)
            shift
            ncpus=$1
            ;;
    -mem | --memory)
            shift
            nram=$1
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



################################################################################

cd ${PROJECT}

# Check and create directories
mkdir -p ${PROJECT}/output
mkdir -p ${PROJECT}/checkpoints
mkdir -p ${PROJECT}/scripts

# Define remaining variables
nFiles=$(ls ${INPUT}/*_R1* | wc -l)

# Create a job script for submission
cat > ${PROJECT}/scripts/alignment_job.sh << EOL
#!/bin/bash
#SBATCH --time=${TIME}  
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1              
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=${ncpus}
#SBATCH --mem=${nram}G               
#SBATCH --job-name=alignment
#SBATCH --output=o.%x_%A_%a.txt    
#SBATCH --error=e.%x_%A_%a.txt      
#SBATCH --mail-user=${EMAIL}
#SBATCH --mail-type=ALL
#SBATCH --array=1-${nFiles}    
###############################################################################
DIR=${INPUT}
mapfile -t FILES < <(find \${DIR} -name "*_R1*")
index=\$((\$SLURM_ARRAY_TASK_ID-1))
inputFile=\${FILES[\$index]}
fileName=\$(basename \$inputFile | awk -F "_R1" '{print \$1}')

bash ${SAMSA2_DIR}/scripts/align.sh -p ${PROJECT} -d \${DIR} -i \${inputFile} -db ${diamond_DB} -a ${aggregation_DB} -c ${common_script} &>> ${PROJECT}/output/step2_\${fileName}.log

EOL

# # Submit the job script
sbatch ${PROJECT}/scripts/alignment_job.sh

