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

#!/bin/bash

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
done
