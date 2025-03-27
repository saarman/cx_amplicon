#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=16                      # Number of CPU cores per node
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

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                     # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                     # Directory to store per-sample VCF files
out_dir="./../cx_amplicon_consensus"               # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/Rep_Genera_Mito.fasta" # Reference FASTA file
region_file="ref_list_ace2.txt"                    # Text file with one region: e.g., ON563187.1:1-712
gene="ace2"                                        # Target gene name for file labeling

# Create output directories if they don't exist
mkdir -p "$out_dir"
mkdir -p "$vcf_dir"

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*

# Read the region of interest from the region file
region=$(cat "$region_file")

# Move into the BAM directory to simplify file handling
cd "$bam_dir"

# Process each BAM file individually
for bam in *.bam; do
  sample=$(basename "$bam" .bam)  # Extract sample name by removing .bam extension

  echo "Processing $sample..."

  # Generate consensus sequence for the specific region
  # --missing N masks low-coverage positions (./. in VCF) with Ns
  samtools faidx "$ref" "$region" | \
  bcftools consensus --missing N "${vcf_dir}/${sample}.vcf" | \
  awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${sample}_${gene}_consensus.fa"
done

# Combine all per-sample consensus FASTAs into one multi-FASTA file
cat "$out_dir"/*_"$gene"_consensus.fa > "$out_dir/all_${gene}_consensus.fa"

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*
