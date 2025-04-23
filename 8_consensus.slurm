#!/bin/bash

#SBATCH --partition=saarman-shared-np
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00
#SBATCH --mem=24576
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name="consensus"

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
vcf_dir="../cx_amplicon_vcf"                      # Directory containing per-sample VCF files
out_dir="../cx_amplicon_consensus"                # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"  # Reference FASTA file

# Step 1: Generate depth-filtered VCFs once
cd "$vcf_dir"
for vcf in *.vcf.gz; do
  sample=$(basename "$vcf" .vcf.gz)
  if [ ! -f "${sample}_DP5.vcf.gz" ]; then
    echo "Filtering $sample to create ${sample}_DP5.vcf.gz"
    bcftools view -Ou "$vcf" | \
      bcftools +setGT -Ou -- -t q -i 'FORMAT/DP<5' -n . | \
      bcftools view -i 'FMT/DP>=2 || GT="./."' -Ov -o ${sample}_DP5.vcf
  bgzip -f ${sample}_DP5.vcf
  tabix -p vcf ${sample}_DP5.vcf.gz
  fi
done

# Step 2: Loop over each gene
for gene in cqm1 ace2 COi; do
  echo "Processing gene: $gene"
  region_file="../cx_amplicon_scripts/ref_list_${gene}.txt"
  region=$(cat "$region_file")
  mkdir -p "${out_dir}/${gene}"
  mkdir -p "${vcf_dir}/${gene}"

  # Loop over samples
  for vcf in *_DP5.vcf.gz; do    
    sample=$(basename "$vcf" _DP5.vcf.gz)

    # From original (unfiltered) VCF (DP2)
    echo "Generating consensus for $sample in $gene for DP2..."
    bcftools view -r "$region" "${sample}.vcf.gz" -Ov -o "${vcf_dir}/${gene}/${sample}_${gene}_DP2.vcf"
    bgzip -f "${vcf_dir}/${gene}/${sample}_${gene}_DP2.vcf"
    tabix -p vcf "${vcf_dir}/${gene}/${sample}_${gene}_DP2.vcf.gz"
    rm -f "${vcf_dir}/${gene}/${sample}_${gene}_DP2.vcf"             # remove vcf to save space
    bcftools consensus --missing N --absent N -f "$ref" "${vcf_dir}/${gene}/${sample}_${gene}_DP2.vcf.gz" | \
      awk -v sample="$sample" -v gene="$gene" '/^>/{print ">"sample"|"gene} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_DP2_consensus.fa"

    # From DP5-filtered VCF
    echo "Generating consensus for $sample in $gene for DP5..."
    bcftools view -r "$region" "${sample}_DP5.vcf.gz" -Ov -o "${vcf_dir}/${gene}/${sample}_${gene}_DP5.vcf"
    bgzip -f "${vcf_dir}/${gene}/${sample}_${gene}_DP5.vcf"
    tabix -p vcf "${vcf_dir}/${gene}/${sample}_${gene}_DP5.vcf.gz"
    rm -f "${vcf_dir}/${gene}/${sample}_${gene}_DP5.vcf"             # remove vcf to save space
    bcftools consensus --missing N --absent N -f "$ref" "${vcf_dir}/${gene}/${sample}_${gene}_DP5.vcf.gz" | \
      awk -v sample="$sample" -v gene="$gene" '/^>/{print ">"sample"|"gene} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_DP5_consensus.fa"
  done

  # Combine all consensus outputs for this gene
  cat ${out_dir}/${gene}/*DP2_consensus.fa > ${out_dir}/all_${gene}_DP2_consensus.fa
  cat ${out_dir}/${gene}/*DP5_consensus.fa > ${out_dir}/all_${gene}_DP5_consensus.fa
done

# Fix permissions
chmod -R g+w ../*
