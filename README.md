# SAMSA2 Pipeline Wrapper
A wrapper for the SAMSA2 pipeline that works with a Slurm job scheduler. This wrapper was designed specifically to work for UBC ARC's Sockeye Cluster

## Installation
To run the samsa2 pipeline clone this directory or download the zip file of this git directory(available on Notion or somewhere).
```
git clone https://github.com/jerryyhe/samsa2.git
```
  
Upon download and extraction of this repository, navigate into the samsa2 directory and set up your Conda environment.
```
conda env create -f samsa2.yml
```

Following creation of your Conda environment, run the additional setup script. This setup script modifies the SAMSA2 conda environment by setting the `PATH` global variable to include the directory containing the samsa2 pipeline when the environment is active and removes the path when enviornment is inactive. Note that if you move the location of the repository that you'll need to rerun this bash script. 
```
bash -i setup.sh
```
***Note that the `-i` flag was used here to run the script in your interactive shell. I find that this is necessary for Sockeye, but I'm not sure why.***  

Activate your Conda environment and you're now ready to run the samsa2 pipeline.
```
conda activate SAMSA2
```

## Usage
This samsa2 wrapper pipeline has two main commands, preprocess and align. Preprocess will trim, filter, merge, and ribodeplete raw sequencing reads. Align will use DIAMOND to align reads to a user-defined DIAMOND database and generate read counts.   

All commands in this samsa2 pipeline require three parameters. First is a project directory specified by `--project_directory` which is the path to the directory where your analysis will take place and is defaulted to the current working directory. In ARC Sockeye, this directory **must** be situated in `scratch` . Additionally, samsa2 also required that path to the directory containing raw sequencing reads specified by `input_directory`. It is highly recommended to have this input situated your project analysis directory in `scratch`. Finally, samsa2 also requires your UBC email for slurm job submission specified by `--email`.  

**Note that all fastq files used as input require `_R1` and `_R2` for paired-end reads and `_R1` for single-end reads.**

#### samsa2.sh
```
Usage: 
    samsa2.sh {preprocess|align}

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
```

#### preprocess
Running `samsa2.sh preprocess` will automatically generate and submit a slurm job script to trim, filter, merge, and ribodeplete raw fastq sequencing files. In addition to project and input file directories, this command also requires a `--min_length` parameter which specifies the minimum read length used for read trimming by Trimmomatic.


#### align
Running `samsa2.sh align` will automatically generate and submit a slurm job script that will perform read alignment and generate read counts. Note that you must have created a DIAMOND database with `diamond makedb` to use as the database parameter for this command. 

