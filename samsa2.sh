#!/bin/bash

export SAMSA2_DIR=$(dirname "$(realpath "$0")")

function usage() {
    cat << USAGE

    Usage: 
        $0 {preprocess|align}

    preprocess:
        -p, --project_directory:            Path of the directory in which you will perform this analysis. Default is the current working directory.
        -i, --input_directory               Input directory. Path to the location of your input reads. Ensure that this directory ONLY contains your input reads and no other files.
        -e, --email                         Your UBC email.
        -m, --min_length                    The minimum read length used for TRIMMOMATIC
        -c, --path_common                   Path to the common.sh script needed to run the pipeline. Default is _____
        -t, --runtime                       The amount of time requested for read preprocessing in hh:mm:ss format. Default is 24:00:00
        -A, --allocation                    The allocation code of your sockeye allocation owner. Default is st-ctropini-1
        -ncpus                              The number of CPUs desired for read preprocessing. Default is 8
        -mem, --memory                      The amount of RAM desired for read preprocessing. Default is 48.

    align:
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

subcommand=$1
shift

case "$subcommand" in
    preprocess)
        echo "Running preprocess_submit.sh script..."
        ${SAMSA2_DIR}/samsa2-master/bash_scripts/preprocess_submit.sh "$@"
        ;;
    subcommand2)
        echo "Running align_submit.sh script..."
        ${SAMSA2_DIR}/align_submit.sh "$@"
        ;;
    *)
        echo "Unknown subcommand: $subcommand"
        usage
        ;;
esac
