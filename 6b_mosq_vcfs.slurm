#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="consensus"                    # Job name for SLURM queue

# Optional email notifications (uncomment if needed)
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load required software modules
module load bwa/2020_03_19
module load samtools/1.16
module load bcftools/1.16
module load htslib

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                     # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                     # Directory to store per-sample VCF files
ref="../cx_amplicon_bwa/ref/Rep_Genera_Mito.fasta" # Reference FASTA file

# Ensure the output directory exists
mkdir -p "$vcf_dir"

# Index the reference FASTA if not already indexed (creates .fai file)
samtools faidx "$ref"

# Apply group write permissions to all directories above (for collaboration or shared environments)
chmod -R g+w ../*

# Move into the BAM directory to simplify file handling
cd "$bam_dir"

# Loop over all BAM files in the directory
for bam in *.bam; do
  sample=$(basename "$bam" .bam)  # Get sample name by stripping the .bam extension

  # Index the BAM file if the index doesn't already exist
  if [ ! -f "${bam}.bai" ] && [ ! -f "${sample}.bai" ]; then
    samtools index "$bam"
  fi

  # Step 1: Generate VCF from BAM using bcftools mpileup filter for mapping quality
  bcftools mpileup -d 2000 -f "$ref" -q 30 "$bam" | \           # MQ filter of 30 (scale is to 60)
  bcftools call -c -Oz -o "${vcf_dir}/${sample}_temp.vcf.gz"    # Output temporary compressed VCF
  tabix -p vcf  "${vcf_dir}/${sample}_temp.vcf.gz"              # Index the compressed VCF

  # Step 2: Filter for regions and depth (DP >= 10) and COI regions
  bcftools view -i 'DP>=10' \
  -r ON563187.1:1-712,KY929304.1:1-1848 "${vcf_dir}/${sample}_temp.vcf.gz" -Oz -o "${vcf_dir}/${sample}_ace2_cqm1.vcf.gz"
  tabix -p vcf "${vcf_dir}/${sample}_ace2_cqm1.vcf.gz"                    # Index the final filtered VCF with tabix

  # Step 3: Cleanup temporary files
  rm "${vcf_dir}/${sample}_temp.vcf.gz" "${vcf_dir}/${sample}_temp.vcf.gz.tbi"

done

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
