#!/bin/bash
#SBATCH --job-name=best_alignments
#SBATCH --output=best_alignments_%j.out
#SBATCH --error=best_alignments_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100  # Adjust based on number of BAM files

module load samtools  

# Set paths
BAM_DIR="/path/to/bam/files"         # Directory containing BAM files
PRIMER_FILE="/path/to/primers.txt"   # File containing primer sequences
OUTPUT_DIR="/path/to/output"         # Output directory for final BAM files

# Get list of BAM files
BAM_FILES=($(ls "$BAM_DIR"/*.bam))
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID - 1]}  # Select BAM file for this job
BASENAME=$(basename "$BAM_FILE" .bam)

# Step 1: Filter out primers and low-quality reads
MAPQ_THRESHOLD=30  # Adjust if needed
samtools view -h "$BAM_FILE" | grep -v -f "$PRIMER_FILE" | samtools view -bq $MAPQ_THRESHOLD -o "$OUTPUT_DIR/${BASENAME}_filtered_q${MAPQ_THRESHOLD}.bam"

# Step 2: Sort BAM by MAPQ (descending) and alignment score
samtools sort -O BAM -t "MAPQ" -o "$OUTPUT_DIR/${BASENAME}_sorted.bam" "$OUTPUT_DIR/${BASENAME}_filtered_q${MAPQ_THRESHOLD}.bam"

# Step 3: Extract the top 100 BEST alignments based on MAPQ
samtools view -h "$OUTPUT_DIR/${BASENAME}_sorted.bam" | awk 'NR<=101 || /^@/' | samtools view -b -o "$OUTPUT_DIR/${BASENAME}_top100.bam"

# Log progress
echo "Processed $BAM_FILE -> $OUTPUT_DIR/${BASENAME}_top100.bam"
