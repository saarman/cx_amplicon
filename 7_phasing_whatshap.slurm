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

# Load modules
module load whatshap/2.6
module load samtools
module load bcftools

# Set directory paths
bam_dir="../cx_amplicon_bwa"
vcf_dir="../cx_amplicon_vcf"
phased_dir="../cx_amplicon_phased"
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"

# Make sure phased directory exists
mkdir -p "$phased_dir"

# Loop through each VCF file
for vcf_file in ${vcf_dir}/*.vcf.gz; do
    sample=$(basename "$vcf_file" .vcf.gz)
    bam_file="${bam_dir}/${sample}.bam"
  #  rg_bam="${bam_dir}/${sample}.rg.bam"

    echo "Phasing $sample"

    # Add read group if not already added
#    if ! samtools view -H "$bam_file" | grep -q "^@RG"; then
#        echo "  No read group found — adding one"
 #       samtools addreplacerg \
 #           -r "@RG\tID:grp1\tSM:${sample}" \
 #           -o "$rg_bam" \
 #           "$bam_file"
 #       samtools index "$rg_bam"
 #   else
  #      rg_bam="$bam_file"
   #     echo "  Read group already present"
  #  fi

    # Run WhatsHap
    whatshap phase \
        -o "${phased_dir}/${sample}.phased.vcf.gz" \
        --reference "$ref" \
        --ignore-read-groups
        "$vcf_file" \
        "$bam_file"
done
