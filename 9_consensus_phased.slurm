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

# Load required software modules, bcftools/1.21 allow us to use --iupac-codes 
module load gcc/11.2.0
module load gcc/11.2.0-cpu
module load gcc/11.2.0-gpu
module load samtools/1.16
module load bcftools/1.21   
module load htslib

# Define input/output paths and filenames                   
out_dir="../cx_amplicon_consensus"                # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"  # Reference FASTA file
vcf_dir="../cx_amplicon_phased"                # Directory containing per-sample phased VCF files

  # Fix permissions
  chmod -R g+w ../*

# Step 1: Generate depth-filtered VCFs once
cd "$vcf_dir"
for vcf in $(bcftools query -l *.phased.vcf.gz | sort | uniq); do
  sample=$(basename "$vcf" .bam)
  if [ ! -f "${sample}_DP5.phased.vcf.gz" ]; then
    echo "Filtering $sample to create ${sample}_DP5.phased.vcf.gz"
    bcftools view -Ou "${sample}.phased.vcf.gz" | \
      bcftools +setGT -Ou -- -t q -i 'FORMAT/DP<5' -n . | \
      bcftools view -i 'FMT/DP>=2 || GT="./."' -Ov -o ${sample}_DP5.phased.vcf
  bgzip -f ${sample}_DP5.phased.vcf
  tabix -p vcf ${sample}_DP5.phased.vcf.gz
  fi
done

# Step 2: Loop over each gene
for gene in cqm1 ace2 COi; do
  echo "Processing phased samples for gene: $gene"
  region_file="../cx_amplicon_scripts/ref_list_${gene}.txt"
  region=$(cat "$region_file")
  mkdir -p "${out_dir}/${gene}"
  mkdir -p "${vcf_dir}/${gene}"

  # Fix permissions
  chmod -R g+w ../*
  
  # Loop over filtered, phased vcf files
  for vcf in *_DP5.phased.vcf.gz; do
    sample=$(basename "$vcf" _DP5.phased.vcf.gz)

    echo "Generating phased consensus for $sample in $gene..."
    bcftools view -r "$region" "$vcf" -Ov -o "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf"
    bgzip -f "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf"
    tabix -f -p vcf "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf.gz"
    rm -f "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf"

    name=$(bcftools query -l "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf.gz")

    # Check for phased GT
    phased=$(bcftools query -f '[%GT]\n' "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf.gz" | grep "|" | wc -l)

    if [ "$phased" -gt 0 ]; then
      echo "Phased genotypes found for $sample. Creating hapA and hapB."

      samtools faidx "$ref" "$region" | \
      bcftools consensus --missing N --absent "-" --sample "$name" --haplotype 1 "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"|hapA"} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_DP5_phased_consensus_hapA.fa"

      samtools faidx "$ref" "$region" | \
      bcftools consensus --missing N --absent "-" --sample "$name" --haplotype 2 "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"|hapB"} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_DP5_phased_consensus_hapB.fa"

    else
      echo "No phased genotypes for $sample. Creating single consensus."
      samtools faidx "$ref" "$region" | \
      bcftools consensus --missing N --absent "-" --iupac-codes --sample "$name" "${vcf_dir}/${gene}/${sample}_${gene}_DP5.phased.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_DP5_phased_consensus.fa"
    fi
  done

  # Combine phased outputs
  cat ${out_dir}/${gene}/*_phased_consensus*.fa > ${out_dir}/all_${gene}_phased_consensus.fa
done

# Fix permissions
chmod -R g+w ../*
