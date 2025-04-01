#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
# #SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="bwa-index"
# #SBATCH --mail-user=emily.calhoun@usu.edu    # email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# Load modules
module load bwa/2020_03_19
#module load bwa-mem/2.2.1

# Replace these variables with the actual paths and filenames
input_ref="ace2_cqm1_coi2.fasta"
output_folder="./../cx_amplicon_bwa/ref/"

# Create the output folder if it doesn't exist, copy fasta file  into output
mkdir -p ${output_folder}
cp $input_ref ${output_folder}/${input_ref}
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon*

# Index reference
cd $output_folder
bwa index -a bwtsw $input_ref    
#bwa-mem2 index Cpip29_nodeb_gfill_mito.fasta

# Permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon*
