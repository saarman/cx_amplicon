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

  # Step 1: Generate VCF from BAM using bcftools mpileup filter for high mapping quality
  bcftools mpileup -d 2000 -f "$ref" -q 50 --ploidy 1 "$bam" | \           # High MQ filter of 50 (scale is to 60)
  bcftools call -c -Oz -o "${vcf_dir}/${sample}_temp.vcf.gz"               # Output temporary compressed VCF
  tabix -p vcf  "${vcf_dir}/${sample}_temp.vcf.gz"                         # Index the compressed VCF

  # Step 2: Filter for regions and depth (DP >= 10) and COI regions
  bcftools view -i 'DP>=10' \
  -r NC_003408.1:1-16697,NC_006853.1:1-16338,NC_008092.1:1-16729,NC_011818.1:1-18266,NC_012920.1:1-16569,NC_013978.1:1-17229,NC_014574.1:1-15587,NC_015079.1:1-14856,NC_020648.1:1-16538,NC_020679.1:1-16352,NC_024867.1:1-16765,NC_024872.1:1-16669,NC_025610.1:1-16822,NC_025611.1:1-16802,NC_025616.1:1-16811,NC_029360.1:1-16793,NC_031863.1:1-17132,NC_037513.1:1-17160,NC_041664.1:1-16804,NC_053523.1:1-16784,MW574355.1:1-16707,NC_072703.1:1-16945,NC_010089.1:1-17446 \
  "${vcf_dir}/${sample}_temp.vcf.gz" -Oz -o "${vcf_dir}/${sample}_coi.vcf.gz"
  tabix -p vcf "${vcf_dir}/${sample}_coi.vcf.gz"                    # Index the final filtered VCF with tabix

  # Step 3: Cleanup temporary files
  rm "${vcf_dir}/${sample}_temp.vcf.gz" "${vcf_dir}/${sample}_temp.vcf.gz.tbi"
done

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
