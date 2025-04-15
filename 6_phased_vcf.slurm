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
min_depth=10                                       # Minimum depth for consensus

# Ensure the output directory exists
mkdir -p "$vcf_dir"
mkdir -p "$phased_dir"

# Index the reference FASTA if not already indexed (creates .fai file)
if [ ! -f "${ref}.fai" ]; then
  samtools faidx "$ref"
fi

# Apply group write permissions to all directories above (for collaboration or shared environments)
chmod -R g+w ../*

# Move into the BAM directory to simplify file handling
cd "$bam_dir"

# Loop over all BAM files in the directory
for bam in *.bam; do
  sample=$(basename "$bam" .bam)  # Get sample name by stripping the .bam extension

  # Always reindex the BAM file to ensure accuracy
  echo "Reindexing BAM file: $bam"
  samtools index "$bam" || { echo "Error indexing BAM file: $bam"; continue; }

  echo "Processing sample: $sample"

  # Step 1: Generate the VCF file with quality filtering using bcftools mpileup
  bcftools mpileup -Ou -f "$ref" -q 30 -Q 30 -d 10000 "$bam" | \
  bcftools call -mv -Oz -o "${vcf_dir}/${sample}.vcf.gz"

  # Step 2: Filter the VCF by minimum per-sample depth
  bcftools view -i "DP>=$min_depth" "${vcf_dir}/${sample}.vcf.gz" -Oz -o "${vcf_dir}/${sample}_filtered.vcf.gz"
  
  # Step 3: Normalize indels to ensure proper representation
  bcftools norm -f "$ref" "${vcf_dir}/${sample}_filtered.vcf.gz" -Oz -o "${vcf_dir}/${sample}_norm.vcf.gz"

  # Step 4: Index the normalized VCF file
  tabix -p vcf "${vcf_dir}/${sample}_norm.vcf.gz"

  # Step 5: Call variants using FreeBayes (haplotype-aware)
  echo "Calling variants for: $sample using FreeBayes"
  freebayes -f "$ref" \
      --ploidy 2 \
      --min-base-quality 30 \
      --min-mapping-quality 30 \
      --min-coverage $min_depth \
      --genotype-qualities \
      "$bam" > "${phased_dir}/${sample}_phased.vcf" || { echo "Error calling variants with FreeBayes for $sample"; continue; }

    # Check if the raw VCF file was created successfully
    if [[ ! -s "${phased_dir}/${sample}_phased.vcf" ]]; then
        echo "Error: Raw VCF file is empty for $sample"
        continue
    fi

    # Compress and index the phased VCF file
    echo "Compressing and indexing ace2_raw VCF: ${sample}_phased.vcf"
    bgzip -c "${phased_dir}/${sample}_phased.vcf" > "${phased_dir}/${sample}_phased.vcf.gz" || { echo "Error compressing phased VCF for $sample"; continue; }
    tabix -p vcf "${phased_dir}/${sample}_phased.vcf.gz" || { echo "Error indexing phased VCF for $sample"; continue; }

done

echo "All samples processed."

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
