#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="cqm1"                         # Job name for SLURM queue

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
phased_dir="../cx_amplicon_phased"                      # Directory to store per-sample phased VCF files
out_dir="../cx_amplicon_consensus"                      # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"        # Reference FASTA file

# Create output directories if they don't exist
mkdir -p "$out_dir"
mkdir -p "${out_dir}/${gene}"

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*

# Move into the VCF directory to simplify file handling
cd "$vcf_dir"

for gene in cqm1 ace2 COi; do
  
  region_file="ref_list_${gene}.txt"
  region=$(cat "$region_file")

  for vcf in ${vcf_dir}/*.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)

    echo "Processing $sample for $gene..."

    bcftools +setGT -Ou -- -t q -i 'FORMAT/DP<5' -n . "$vcf" | \
    bcftools view -i 'FMT/DP>=2 || GT="./."' -r "$region" -Ou | \
    bcftools consensus --missing N -f "$ref" | \
    awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_consensus.fa"

  done
done

# Combine all per-sample consensus FASTAs into one multi-FASTA file
cat ${out_dir}/${gene}/*_consensus.fa > ${out_dir}/all_${gene}_consensus.fa

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*
