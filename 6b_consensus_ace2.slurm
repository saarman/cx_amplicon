#!/bin/bash

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576       # memory given in MB
#SBATCH --nodes=1         # number of nodes
#SBATCH --ntasks-per-node=16
#SBATCH --job-name="consensus"

# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load modules
module load bwa/2020_03_19
module load samtools/1.16
module load bcftools/1.16

# Paths and filenames
bam_dir="./../cx_amplicon_bwa"
vcf_dir="./../cx_amplicon_vcf"
out_dir="./../cx_amplicon_consensus"
ref="../cx_amplicon_bwa/ref/Rep_Genera_Mito.fasta"
region_file="ref_list_ace2.txt"  # File contains one line: ON563187.1:1-712
gene="ace2"

# Make output directory if it doesn't exist
mkdir -p "$out_dir"
chmod -R g+w ../*  # Permissions

# Read region from file
region=$(cat "$region_file")

# Loop through each .bam file in the bam_dir
cd "$bam_dir"
for bam in *.bam; do
  sample=$(basename "$bam" .bam)

  # Ensure BAM is indexed
  if [ ! -f "${bam}.bai" ] && [ ! -f "${sample}.bai" ]; then
    samtools index "$bam"
  fi

  echo "Processing $sample..."

  # Create per-sample VCF
  bcftools mpileup -d 1000 -f "$ref" "$bam" | \
  bcftools call -c -Ov -o "${vcf_dir}/${sample}_${gene}.vcf"

  # Compress and index the VCF
  bgzip -f "${vcf_dir}/${sample}_${gene}.vcf"
  tabix -p vcf "${vcf_dir}/${sample}_${gene}.vcf.gz"

  # Generate consensus for the region using samtools faidx
  samtools faidx "$ref" "$region" | \
  bcftools consensus "${vcf_dir}/${sample}_${gene}.vcf.gz" | \
  awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${sample}_${gene}_consensus.fa"
done

# Concatenate all sample consensus FASTAs into one file
cat "$out_dir"/*_"$gene"_consensus.fa > "$out_dir/all_${gene}_consensus.fa"
chmod -R g+w ../*  # Permissions
