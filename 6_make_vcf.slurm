#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="vcf"                          # Job name for SLURM queue

# Optional email notifications (uncomment if needed)
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load required software modules
module load freebayes/1.3.4
module load bwa/2020_03_19
module load samtools/1.16
module load bcftools/1.16
module load htslib

# Working directory /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_scripts

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                     # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                     # Directory to store per-sample VCF files
phased_dir="./../cx_amplicon_phased"                  # Directory to store per-sample phased VCF files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"   # Reference FASTA file

# Ensure the output directory exists
mkdir -p "$vcf_dir"

# Index the reference FASTA (creates .fai file)
samtools faidx "$ref"

# Apply group write permissions to all directories above (for collaboration or shared environments)
chmod -R g+w ../*

# Move into the BAM directory to simplify file handling
cd "$bam_dir"

# Loop over all BAM files in the directory
for bam in *.bam; do
  sample=$(basename "$bam" .bam)  # Get sample name by stripping the .bam extension

  # Index the BAM file
  samtools index "$bam"

  echo "Processing sample: $sample"

  # Step 1: Generate the VCF file with quality filtering using bcftools mpileup
  # includes sites with zero coverage, marked with no genotype (./.) and missing depth
  bcftools mpileup -Ou -f "$ref" -q 30 -Q 30 -d 10000 -a FORMAT/DP "$bam" | \
  bcftools call -m --ploidy 2 -Ou | \
  bcftools +setGT -Ou -- -t q -i 'FORMAT/DP<2' -n . | \
  bcftools norm -f "$ref" -Ou | \
  bcftools view -i 'DP>=1 || GT="./."' -Oz -o "${vcf_dir}/${sample}.vcf.gz"

  # Step 2: Index the VCF file, forcing overwrite
  tabix -p vcf "${vcf_dir}/${sample}.vcf.gz" -f

  # CHECK: Are we creating correct (./.) with no data? Use this...
  # bcftools query -f '%GT\n' ../cx_amplicon_vcf/*_L001.vcf.gz | grep '\./\.'  
  # bcftools query -f '%GT\n' ../cx_amplicon_vcf/*_norm.vcf.gz | grep '\./\.'
  
  # Future Improvement: Call variants using FreeBayes (haplotype-aware)
  # SKIP (FOR NOW), THE REPORT ALL HAPLOTYPES IS NOT WORKING
  # MAYBE CAN USE FREEBAYES AFTER THE VCF FILE IS ALREADY CREATED?
  # MAYBE SHOULD USE A DIFFERENT PHASING SOFTWARE
  #echo "Calling variants for: $sample using FreeBayes"
  #freebayes -f "$ref" \
  # --ploidy 2 \
  # --min-base-quality 30 \
  # --min-mapping-quality 30 \
  # --min-coverage $min_depth \
  # --genotype-qualities \
  # --report-monomorphic \
  # --report-all-haplotype-alleles \
  # "$bam" > "${phased_dir}/${sample}_phased.vcf" || { echo "Error calling variants with FreeBayes for $sample"; continue; }
  #
  # Check if the raw VCF file was created successfully
  #if [[ ! -s "${phased_dir}/${sample}_phased.vcf" ]]; then
  #   echo "Error: Raw VCF file is empty for $sample"
  #   continue
  #fi
  #
  # Compress and index the phased VCF file
  #bgzip -c "${phased_dir}/${sample}_phased.vcf" > "${phased_dir}/${sample}_phased.vcf.gz" || { echo "Error compressing phased VCF for $sample"; continue; }
  #tabix -p vcf "${phased_dir}/${sample}_phased.vcf.gz" || { echo "Error indexing phased VCF for $sample"; continue; }
done

echo "All samples processed."

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
