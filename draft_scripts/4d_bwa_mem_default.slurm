#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=bwa-mem
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL
# #SBATCH --mail-user=emily.calhoun@usu.edu

# Load modules
module load bwa/2020_03_19
module load samtools/1.16

# Run the Perl script with the input files
perl ./4d_bwa_mem_default.pl ../../cx_amplicon_raw/*R1_001.fastq.gz  

# Permissions
chmod -R g+w ../*  

# Resource
# This is calling a perl script that uses ForkManager
