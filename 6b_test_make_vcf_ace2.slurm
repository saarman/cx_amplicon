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

# NOTE: Working directory /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_scripts

# Load required software modules
module load bwa/2020_03_19
module load samtools/1.16
module load bcftools/1.16
module load htslib

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                      # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                       # Directory to store per-sample VCF files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"     # Reference FASTA file

# Ensure the output directory exists
mkdir -p "$vcf_dir"

# Index the reference FASTA file each time (to account for updates)
echo "Indexing reference FASTA: $ref"
samtools faidx "$ref" || { echo "Error indexing reference"; exit 1; }

# Apply group write permissions to all directories above (for collaboration or shared environments)
chmod -R g+w ../*

# Move into the BAM directory to simplify file handling
cd "$bam_dir" || { echo "Error: Cannot change to BAM directory: $bam_dir"; exit 1; }

# Loop over all BAM files matching the pattern
for bam in B053-UT-M07101-240702*.bam B054-UT-M07101-240702*.bam B292-UT-M70330-240718*.bam B002-UT-M07101-240702*.bam B002-UT-M70330-240718*.bam; do
    # Check if the BAM file exists
    if [[ ! -f "$bam" ]]; then
        echo "Warning: No BAM files matching pattern: $bam"
        continue
    fi

    sample=$(basename "$bam" .bam)  # Get sample name by stripping the .bam extension
    echo "Processing sample: $sample"

    # Always reindex the BAM file to ensure accuracy
    echo "Reindexing BAM file: $bam"
    samtools index "$bam" || { echo "Error indexing BAM file: $bam"; continue; }

    # 1. Generate genotype likelihoods (BCF format is recommended for phasing)
    echo "Generating genotype likelihoods for: $sample"
    bcftools mpileup -d 2000 -f "$ref" -q 20 -Q 13 "$bam" -Ou | \
    bcftools call -m -v -Oz -o "${vcf_dir}/${sample}_raw.vcf.gz" || { echo "Error in mpileup or call for $sample"; continue; }

    # 2. Index the raw VCF file
    echo "Indexing raw VCF: ${sample}_raw.vcf.gz"
    bcftools index "${vcf_dir}/${sample}_raw.vcf.gz" || { echo "Error indexing raw VCF for $sample"; continue; }

    # 3. Phase the haplotypes
    echo "Phasing haplotypes for: $sample"
    bcftools phase -O z -o "${vcf_dir}/${sample}_phased.vcf.gz" "${vcf_dir}/${sample}_raw.vcf.gz" || { echo "Error phasing for $sample"; continue; }

    # 4. Filter for depth (DP >= 10)
    echo "Filtering phased VCF for depth (DP >= 10): $sample"
    bcftools view -i 'DP>=10' "${vcf_dir}/${sample}_phased.vcf.gz" -Oz -o "${vcf_dir}/${sample}_filtered_phased.vcf.gz" || { echo "Error filtering VCF for $sample"; continue; }

    # 5. Index the filtered phased VCF file
    echo "Indexing filtered phased VCF: ${sample}_filtered_phased.vcf.gz"
    tabix -p vcf "${vcf_dir}/${sample}_filtered_phased.vcf.gz" || { echo "Error indexing filtered VCF for $sample"; continue; }

    # 6. Generate a consensus sequence from the phased VCF
    echo "Generating consensus sequence for: $sample"
    bcftools consensus -f "$ref" "${vcf_dir}/${sample}_filtered_phased.vcf.gz" > "${vcf_dir}/${sample}_phased_consensus.fa" || { echo "Error generating consensus for $sample"; continue; }

    # Cleanup intermediate VCF files to save space
    echo "Cleaning up intermediate files for: $sample"
    rm -f "${vcf_dir}/${sample}_raw.vcf.gz" "${vcf_dir}/${sample}_raw.vcf.gz.csi" \
          "${vcf_dir}/${sample}_phased.vcf.gz" "${vcf_dir}/${sample}_phased.vcf.gz.csi"
    
    echo "Completed processing for: $sample"
done

echo "All BAM files processed successfully."

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
