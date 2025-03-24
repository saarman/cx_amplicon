#!/bin/bash
#SBATCH --time=6:00:00
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
BAM_DIR="./../cx_amplicon_bwa"                      # Directory containing BAM files
PRIMER_FILE="../cx_amplicon_bwa/ref/primers.fasta"  # File containing primer sequences
OUTPUT_DIR="./../cx_amplicon_30filter"           # Output directory for final BAM files
MAPQ_THRESHOLD=30                                   # Change if needed//checked and 30 looks good
TOP_N=100                                           # Number of best alignments to keep

mkdir -p "$OUTPUT_DIR"

# Loop through all BAM files
for BAM_FILE in "$BAM_DIR"/*.bam; do
    BASENAME=$(basename "$BAM_FILE" .bam)

    echo "Processing $BASENAME..."

    # Step 1: Filter out primers and low-quality alignments
    samtools view -h "$BAM_FILE" | grep -v -f "$PRIMER_FILE" | \
        samtools view -bq "$MAPQ_THRESHOLD" -o "$OUTPUT_DIR/${BASENAME}_filtered_q${MAPQ_THRESHOLD}.bam"

    # Check if the filtered BAM file is empty
    if [[ ! -s "$OUTPUT_DIR/${BASENAME}_filtered_q${MAPQ_THRESHOLD}.bam" ]]; then
        echo "No valid reads after filtering for $BASENAME. Skipping..."
        continue
    fi

    # Step 2: Sort by MAPQ in descending order (using samtools sort)
    samtools sort -n -o "$OUTPUT_DIR/${BASENAME}_sorted.bam" "$OUTPUT_DIR/${BASENAME}_filtered_q${MAPQ_THRESHOLD}.bam"

    echo "Sorted BAM saved as: $OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Step 3: Extract top 100 best alignments
    samtools view -h "$OUTPUT_DIR/${BASENAME}_sorted.bam" | head -n $((TOP_N + 1)) | \
        samtools view -b -o "$OUTPUT_DIR/${BASENAME}_top${TOP_N}.bam"

    echo "Top $TOP_N alignments saved as: $OUTPUT_DIR/${BASENAME}_top${TOP_N}.bam"
done

echo "All BAM processing complete."
