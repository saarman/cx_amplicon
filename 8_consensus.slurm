#!/bin/bash

#SBATCH --partition=saarman-shared-np
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00
#SBATCH --mem=24576
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name="cqm1"

# Optional email notifications
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load required software modules
module load samtools/1.16
module load bcftools/1.16
module load htslib

# Define input/output paths and filenames
vcf_dir="../cx_amplicon_vcf"            # Directory containing per-sample VCF files
out_dir="../cx_amplicon_consensus"       # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi_clean.fasta"  # Reference FASTA file

# Move into the VCF directory
cd "$vcf_dir"

# Filter VCF for depth as you wish
for vcf in *.vcf.gz; do
  sample=$(basename "$vcf" .vcf.gz)
  
  echo "Filtering $sample for depth..."

  # Stream and process, output VCF filtering for depth
  bcftools view -Ou "$vcf" | \
  bcftools +setGT -Ou -- -t q -i 'FORMAT/DP<5' -n . | \
  bcftools view -i 'FMT/DP>=2 || GT="./."' -Ov -o ${sample}_DP5.vcf
  
  # Zip and index
  bgzip -f ${sample}_DP5.vcf
  tabix -p vcf ${sample}_DP5.vcf.gz
done

for gene in cqm1 ace2 COi; do
  # Create output directories if they don't exist
  mkdir -p "${out_dir}/${gene}"

  # Define region
  region_file="../cx_amplicon_scripts/ref_list_${gene}.txt"
  region=$(cat "$region_file")

  for vcf in *.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)

    echo "Processing $sample for $gene..."

    # Extract region and make consensus, choose what input here!!! Filter for depth beyond 2?
    bcftools view -r "$region" ${sample}.vcf.gz -Ov -o temp_region.vcf

    # Compress and index
    bgzip -f temp_region.vcf
    tabix -p vcf temp_region.vcf.gz

    bcftools consensus --missing N --absent N -f "$ref" temp_region.vcf.gz - | \
    awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_consensus.fa"

    # Clean up
    rm temp*.vcf.gz temp*.vcf.gz.tbi
  done

  # Combine all per-sample consensus FASTAs for this gene
  cat ${out_dir}/${gene}/*_consensus.fa > ${out_dir}/all_${gene}_consensus.fa

done

# Ensure group write permissions for all relevant directories
chmod -R g+w ../*
