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
vcf_dir="../cx_amplicon_vcf"                      # Directory containing per-sample VCF files
out_dir="../cx_amplicon_consensus"                # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"  # Reference FASTA file

# Step 1: Generate depth-filtered VCFs (DP3)
cd "$vcf_dir"
for vcf in $(bcftools query -l *.vcf.gz | sort | uniq); do
  sample=$(basename "$vcf" .bam)
  if [ ! -f "${sample}_DP3.vcf.gz" ]; then
    echo "Filtering $sample to create ${sample}_DP3.vcf.gz"

    # Looser DP3: positions with DP < 3 → missing
    bcftools view -Ou "${sample}.vcf.gz" | \
      bcftools +setGT -Ou -- -t q -i 'FORMAT/DP<3' -n . | \
      bcftools view -i 'FMT/DP>=2 || GT="./."' -Ov -o ${sample}_DP3.vcf

    bgzip -f ${sample}_DP3.vcf
    tabix -p vcf ${sample}_DP3.vcf.gz
  fi
done

# Step 2: Loop over each gene for DP3 consensus
for gene in cqm1 ace2 COi; do
  echo "Processing gene: $gene for DP3 consensus"
  region_file="../cx_amplicon_scripts/ref_list_${gene}.txt"
  region=$(cat "$region_file")
  mkdir -p "${out_dir}/${gene}"
  mkdir -p "${vcf_dir}/${gene}"

  # Loop over DP3-filtered VCFs
  for vcf in *_DP3.vcf.gz; do    
    sample=$(basename "$vcf" _DP3.vcf.gz)

    echo "Generating DP3 consensus for $sample in $gene..."
    bcftools view -r "$region" "$vcf" -Ov -o "${vcf_dir}/${gene}/${sample}_${gene}_DP3.vcf"
    bgzip -f "${vcf_dir}/${gene}/${sample}_${gene}_DP3.vcf"
    tabix -f -p vcf "${vcf_dir}/${gene}/${sample}_${gene}_DP3.vcf.gz"
    rm -f "${vcf_dir}/${gene}/${sample}_${gene}_DP3.vcf"  # remove temporary VCF

    name=$(bcftools query -l "${vcf_dir}/${gene}/${sample}_${gene}_DP3.vcf.gz")
    samtools faidx "$ref" "$region" | \
    bcftools consensus --missing N --absent "-" --iupac-codes --sample "$name" "${vcf_dir}/${gene}/${sample}_${gene}_DP3.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_DP3_consensus.fa"
  done

  # Combine all DP3 consensus sequences for this gene
  cat ${out_dir}/${gene}/*_DP3_consensus.fa > ${out_dir}/all_${gene}_DP3_consensus.fa
done

# Fix permissions
chmod -R g+w ../*

echo "DP3 consensus sequences generated for all genes."


