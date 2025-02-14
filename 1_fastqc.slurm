#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="fastqc"
#SBATCH --mail-user=emily.calhoun@usu.edu   # insert email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# Load Module
module load fastqc/0.12.1

# Input and output files
input_folder="/uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_processed/"
output_folder="/uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_fastqc/"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_ddRAD*

# Run FastQC on all FASTQ files in the input folder
for fastq_file in "$input_folder"/*.fq*; do
    # Check if the file is a regular file
    if [ -f "$fastq_file" ]; then
        # Extract the file name without the path and extension
        base_name=$(basename -s .fq "$fastq_file")

        # Run FastQC on the current FASTQ file
        fastqc "$fastq_file" --outdir="$output_folder"

        # Print a message indicating completion for this file
        echo "FastQC completed for $base_name"
    fi
done
