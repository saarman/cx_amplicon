#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="phasing"                          # Job name for SLURM queue

# Optional email notifications (uncomment if needed)
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load required modules
module load whatshap/2.6
module load bcftools
module load samtools

# Define directories
bam_dir="./../cx_amplicon_bwa"
vcf_dir="./../cx_amplicon_vcf"
phased_dir="./../cx_amplicon_phased"
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"

# Make sure output directory exists
mkdir -p "$phased_dir"

# Index reference if needed
samtools faidx "$ref"

# Loop through each sample VCF and perform phasing
cd "$vcf_dir"
for vcf in *.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)
    bam="${bam_dir}/${sample}.bam"
    
    echo "Phasing $sample"
    
    whatshap phase \
        --reference "$ref" \
        --output "${phased_dir}/${sample}_phased.vcf" \
        "$vcf" "$bam"
    
    # Compress and index phased VCF
    bgzip -f "${phased_dir}/${sample}_phased.vcf"
    tabix -p vcf "${phased_dir}/${sample}_phased.vcf.gz"
done

echo "Phasing complete."
