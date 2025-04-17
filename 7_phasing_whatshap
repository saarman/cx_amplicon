#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="vcf"                          # Job name for SLURM queue

# Optional email notifications (uncomment if needed)
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL



# Paths
CONTAINER="/shared/path/to/whatshap_latest.sif"   # Shared Singularity container location
REF="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"
VCF_DIR="../cx_amplicon_vcf"
BAM_DIR="../cx_amplicon_bwa"
PHASED_DIR="../cx_amplicon_phased"
SAMPLES=($(ls ${VCF_DIR}/*.vcf.gz))  # assumes 1 VCF per sample

# Get sample based on array ID
VCF_FILE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$VCF_FILE" .vcf.gz)
BAM_FILE="${BAM_DIR}/${BASENAME}.bam"
OUT_FILE="${PHASED_DIR}/${BASENAME}_phased.vcf"

# Make output directory if needed
mkdir -p "$PHASED_DIR"
mkdir -p logs

# Step 1: Pull Singularity container with whatshap (if not already done)
singularity pull --name whatshap_latest.sif docker://biocontainers/whatshap:1.4.1--pyhdfd78af_0 || {
    echo "Error pulling Singularity container. Exiting."
    exit 1
}

# Step 2: Set up and activate conda environment inside Singularity container
singularity exec --bind $(pwd):/mnt $CONTAINER bash -c "
    # Create conda environment (only if not already created)
    conda create -n whatshap_env -c bioconda -c conda-forge whatshap -y

    # Activate the conda environment
    source activate whatshap_env

    # Run whatshap for phasing
    whatshap phase --reference /mnt/$REF /mnt/$VCF_FILE /mnt/$BAM_FILE --output /mnt/$OUT_FILE --indels
"

# Step 3: Check if the output file exists and move to final directory
if [[ ! -s "$OUT_FILE" ]]; then
    echo "Error: Phased VCF file is empty for $BASENAME"
    exit 1
fi

# Step 4: Log success and compress the phased VCF file
bgzip -c "$OUT_FILE" > "${PHASED_DIR}/${BASENAME}_phased.vcf.gz"
tabix -p vcf "${PHASED_DIR}/${BASENAME}_phased.vcf.gz"

echo "Phasing completed for $BASENAME"
