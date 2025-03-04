#!/bin/sh
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu

# Load modules
module load bwa/2020_03_19
module load samtools/1.16

# Change to the directory where the input data is located
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_processed

# Run the Perl script with the input files
perl /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_scripts/4b_bwa_mem.pl *.fq.gz

# Permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon*

# Resource
# https://catchenlab.life.illinois.edu/stacks/manual/#refmapbyhand
