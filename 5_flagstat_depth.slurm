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

# Define directories
bam_dir="./../cx_amplicon_bwa"
depth_dir="${bam_dir}/depth"

# Create depth output directory if it doesn't exist
mkdir -p "$depth_dir"

# Output spreadsheet file
output_spreadsheet="flagstat_results.tsv"

# Create a header for the TSV file
echo -e "Sample\tTotal_Reads\tMapped_Reads" > "$output_spreadsheet"

# Loop through all BAM files in the directory
for bam_file in "$bam_dir"/*.bam; do
    # Extract the sample name from the file name
    sample_name=$(basename "$bam_file" .bam)

    # Calculate flagstat metrics
    flagstat_output=$(samtools flagstat "$bam_file")

    # Extract relevant metrics
    total_reads=$(echo "$flagstat_output" | grep "(QC-passed reads + QC-failed reads)" | awk '{print $1}')
    mapped_reads=$(echo "$flagstat_output" |grep "+ 0 mapped" | awk '{print $1}')

    # Append the results to the spreadsheet
    echo -e "$sample_name\t$total_reads\t$mapped_reads" >> "$output_spreadsheet"

    # --- depth output ---
    samtools depth -a "$bam_file" > "${depth_dir}/${sample_name}.depth.txt"
done

# After the loop: Compute mean depth per position across all samples
echo "Computing mean depth per position..."

# Create a combined file: chrom, pos, mean_depth
# Assumes all depth files are in same 3-column format: CHROM POS DEPTH

cd "$depth_dir"

# Use paste to align all depths by position
# First extract positions from the first file
first_file=$(ls *.depth.txt | head -n 1)
cut -f1,2 "$first_file" > all_positions.txt

# Add all depth columns to a temporary matrix
for f in *.depth.txt; do
    cut -f3 "$f" > "onlydepth_${f}"
done

# Combine all depths side-by-side
paste all_positions.txt onlydepth_*.depth.txt > depth_matrix.tsv

# Calculate mean depth per position
awk '{
    sum=0; count=0;
    for(i=3; i<=NF; i++) {
        if ($i != "NA" && $i != ".") {
            sum+=$i; count++;
        }
    }
    mean = (count > 0) ? sum/count : "NA";
    print $1, $2, mean
}' depth_matrix.tsv > mean_depth_per_position.tsv

# Clean up
rm onlydepth_*.depth.txt all_positions.txt depth_matrix.tsv

echo "Mean depth written to: $depth_dir/mean_depth_per_position.tsv"

