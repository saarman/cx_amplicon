#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="ace2"                         # Job name for SLURM queue

# Optional email notifications (uncomment if needed)
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load required software modules
module load samtools/1.16
module load bcftools/1.16

# Define input/output paths and filenames
vcf_dir="../cx_amplicon_vcf"                            # Directory to store per-sample VCF files
out_dir="../cx_amplicon_consensus"                      # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"        # Reference FASTA file
region_file="../cx_amplicon_scripts/ref_list_ace2.txt"  # Text file with one region: e.g., ON563187.1:1-712
gene="ace2"                                             # Target gene name for file labeling

# Create output directories if they don't exist
mkdir -p "$out_dir"
mkdir -p "${out_dir}/${gene}"

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*

# Read the region of interest from the region file
region=$(cat "$region_file")

# Move into the VCF directory to simplify file handling
cd "$vcf_dir"

# Process each VCF file individually
for vcf in *_norm.vcf.gz; do
  sample=$(basename "$vcf" _norm.vcf.gz)  # Extract sample name by removing .vcf extension

  echo "Processing $sample..."

  # Generate consensus, --missing N masks low-coverage positions (./. in VCF) with Ns 
  #bcftools view -i 'DP>=10' "${sample}_norm.vcf.gz"        # gave up extra filter, should aready be done in step 6
  #awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_consensus.fa"
  samtools faidx "$ref" "$region" | bcftools consensus --missing N "${sample}_norm.vcf.gz" | \
  awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_consensus.fa"
done

# Combine all per-sample consensus FASTAs into one multi-FASTA file
cat ${out_dir}/${gene}/*.fa > ${out_dir}/all_${gene}_consensus.fa

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*
