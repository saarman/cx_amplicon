#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20         
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=bwa-mem2
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL
# #SBATCH --mail-user=emily.calhoun@usu.edu

module load samtools/1.16  

# Set paths
BAM_DIR="./../cx_amplicon_bwa"       # Directory containing BAM files
PRIMER_FILE="../cx_amplicon_bwa/ref/primers.fasta"  # File containing primer sequences
OUTPUT_DIR="./../cx_amplicon_filtered"         # Output directory for final BAM files

# Loop through all BAM files in BAM_DIR
for BAM_FILE in "$BAM_DIR"/*.bam; do
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
done
