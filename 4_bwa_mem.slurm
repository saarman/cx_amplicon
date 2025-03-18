#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=bwa-mem2
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL
# #SBATCH --mail-user=emily.calhoun@usu.edu

# Load modules
module load bwa/2020_03_19
module load samtools/1.16
module load samclip

# Run the Perl script with the input files
perl /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_scripts/4_bwa_mem.pl ../../cx_amplicon_raw/*R1_001.fastq.gz    #tried changing from _001 to R1_001

# Permissions
# chmod -R g+w ../*  #commented this out because i think this is messing the script up 

# Resource
# This is calling a perl script that uses ForkManager
