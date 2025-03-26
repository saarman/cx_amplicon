#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="bwa"
# #SBATCH --mail-user=emily.calhoun@usu.edu   # email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load modules
module load bwa/2020_03_19
module load samtools/1.16

# Replace these variables with the actual paths and filenames
bam_dir="./../cx_amplicon_bwa"
out_dir="./../cx_amplicon_consensus"
gene="COi"
scaffolds="COi_scaffolds.list"

# make directory if not already there
mkdir out_dir      # Make directory
chmod -R g+w ../*  # Permissions

# loop through each .bam file in the bam_dir
cd ${bam_dir}
for bam in *.bam; do
  sample=$(basename "$bam" .bam)

  # Ensure BAM is indexed
  if [ ! -f "${bam}.bai" ] && [ ! -f "${sample}.bai" ]; then
    samtools index "$bam"
  fi

  echo "Processing $sample..."

  samtools mpileup -uf reference.fasta "$bam" | \
  bcftools call -c | \
  bcftools filter -s LowQual -e '%QUAL<20 || DP<10' | \
  bcftools view -R ${scaffolds} | \
  bcftools consensus -f ${ref} | \
  awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > ${outdir}/${sample}_${gene}_consensus.fa
done

# Concatenate all coi output from all samples into one file
cat ${outdir}/*_${gene}_consensus.fa > all_${gene}_consensus.fa
chmod -R g+w ../*  # Permissions
