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

  # Generate VCF from BAM using bcftools mpileup + call, no filtering
  # -d 1000 sets a high depth cap suitable for amplicon sequencing
  bcftools mpileup -d 2000 -f "$ref" -q 20 "$bam" | \
  bcftools call -c -Ov -o "${vcf_dir}/${sample}_temp.vcf"

  # 2. Filter for depth (DP >= 30)
  bcftools view -i 'DP>=30' "${vcf_dir}/${sample}_temp.vcf" -Ov -o "${vcf_dir}/${sample}.vcf"

  # Compressing and indexing ${sample}.vcf...
  bgzip -f "${vcf_dir}/${sample}.vcf"
  tabix -p vcf "${vcf_dir}/${sample}.vcf.gz"
done

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
