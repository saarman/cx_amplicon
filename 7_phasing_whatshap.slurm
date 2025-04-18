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
module load singularity/4.1.1
module load samtools/1.16

# Set directories
bam_dir="./../cx_amplicon_bwa"                    # Directory with BAM files
vcf_dir="./../cx_amplicon_vcf"                    # Directory with input VCFs
phased_dir="./../cx_amplicon_phased"              # Output dir for phased VCFs
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"  # Reference FASTA
container="./containers/whatshap_1.4.1--pyhdfd78af_0.sif"  # Path to whatshap container

# Check if container exists
if [[ ! -f "$container" ]]; then
    echo "ERROR: Singularity container not found at: $container"
    echo "Please pull it manually with:"
    echo "singularity pull docker://quay.io/biocontainers/whatshap:1.4.1--pyhdfd78af_0"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$phased_dir"

# Loop over each sample VCF
for vcf in "$vcf_dir"/*.vcf.gz; do
  sample=$(basename "$vcf" .vcf.gz)

  echo "Phasing sample: $sample"

  # Run Whatshap phasing
  singularity exec "$container" whatshap phase \
    -o "${phased_dir}/${sample}_phased.vcf.gz" \
    --reference "$ref" \
    "$vcf" \
    "${bam_dir}/${sample}.bam"

  # Index phased VCF
  singularity exec "$container" bcftools index "${phased_dir}/${sample}_phased.vcf.gz"

done

echo "Phasing completed for all samples."
