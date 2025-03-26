#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="consensus"
# #SBATCH --mail-user=emily.calhoun@usu.edu   # email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load modules
module load bwa/2020_03_19
module load samtools/1.16
module load bcftools/1.21

# Paths and filenames
bam_dir="./../cx_amplicon_bwa"
out_dir="./../cx_amplicon_consensus"
ref="../cx_amplicon_bwa/ref/Rep_Genera_Mito.fasta"
scaffolds="ref_list_cqm1.txt"
gene="cqm1"

# Make output directory if it doesn't exist
mkdir -p "$out_dir"
chmod -R g+w ../*  # Permissions

# Loop through each .bam file in the bam_dir
cd "$bam_dir"
for bam in *.bam; do
  sample=$(basename "$bam" .bam)

  # Ensure BAM is indexed
  if [ ! -f "${bam}.bai" ] && [ ! -f "${sample}.bai" ]; then
    samtools index "$bam"
  fi

  echo "Processing $sample..."

  samtools mpileup -uf "$ref" "$bam" | \
  bcftools call -c | \
  bcftools filter -s LowQual -e '%QUAL<20 || DP<10' | \
  bcftools view -R "$scaffolds" | \
  bcftools consensus -f "$ref" | \
  awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${out_dir}/${sample}_${gene}_consensus.fa"
done

# Concatenate all sample consensus FASTAs into one file
cat "$out_dir"/*_"$gene"_consensus.fa > "$out_dir/all_${gene}_consensus.fa"
chmod -R g+w ../*  # Permissions
