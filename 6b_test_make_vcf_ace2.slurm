#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="vcf"                    # Job name for SLURM queue

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

# Working directory /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_scripts

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                     # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                     # Directory to store per-sample VCF files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"  # Reference FASTA file

# Ensure the output directory exists
mkdir -p "$vcf_dir"

# Index the reference FASTA if not already indexed (creates .fai file)
samtools faidx "$ref"

# Apply group write permissions to all directories above (for collaboration or shared environments)
chmod -R g+w ../*

# Move into the BAM directory to simplify file handling
cd "$bam_dir"

# Loop over all BAM files in the directory
for bam in `echo B053-UT-M07101-240702*.bam B054-UT-M07101-240702*.bam B292-UT-M70330-240718*.bam B002-UT-M07101-240702*.bam B002-UT-M70330-240718*bam`; do
  sample=$(basename "$bam" .bam)  # Get sample name by stripping the .bam extension

  # Index the BAM file if the index doesn't already exist
  samtools index "$bam"

  # 1. Generate genotype likelihoods (BCF format is recommended for phasing)
  bcftools mpileup -d 2000 -f "$ref" -q 20 -Q 13 "$bam" -Ou | \
  bcftools call -m -v -Oz -o "${vcf_dir}/${sample}_raw.vcf.gz"

  # 2. Index the raw VCF file
  bcftools index "${vcf_dir}/${sample}_raw.vcf.gz"

  # 3. Phase the haplotypes
  bcftools phase -O z -o "${vcf_dir}/${sample}_phased.vcf.gz" "${vcf_dir}/${sample}_raw.vcf.gz"

  # 4. Filter for depth (DP >= 10)
  bcftools view -i 'DP>=10' "${vcf_dir}/${sample}_phased.vcf.gz" -Oz -o "${vcf_dir}/${sample}_filtered_phased.vcf.gz"

  # 5. Index the filtered phased VCF file
  tabix -p vcf "${vcf_dir}/${sample}_filtered_phased.vcf.gz"

  # 6. Generate a consensus sequence from the phased VCF
  cat "$ref" | bcftools consensus "${vcf_dir}/${sample}_filtered_phased.vcf.gz" > "${vcf_dir}/${sample}_phased_consensus.fa"

  # Cleanup temporary files
  rm "${vcf_dir}/${sample}_temp.vcf.gz" "${vcf_dir}/${sample}_temp.vcf.gz.tbi"
done

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
