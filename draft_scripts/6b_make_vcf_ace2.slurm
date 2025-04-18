#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="freebayes"                          # Job name for SLURM queue

# Optional email notifications (uncomment if needed)
# #SBATCH --mail-user=emily.calhoun@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Working directory /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_scripts

# Load required software modules
module load freebayes/1.3.4
module load samtools/1.16
module load bcftools/1.16
module load htslib

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                      # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                       # Directory to store per-sample VCF files
con_dir="./../cx_amplicon_consensus"                       # Directory to store per-sample consensus
hap_dir="./../cx_amplicon_haplotype"                       # Directory to store per-sample haplotypes
ref="../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta"     # Reference FASTA file

# Ensure the output directory exists
mkdir -p "$vcf_dir"
mkdir -p "$con_dir"
mkdir -p "$hap_dir"

# Index the reference FASTA file each time (to account for updates)
echo "Indexing reference FASTA: $ref"
samtools faidx "$ref" || { echo "Error indexing reference"; exit 1; }

# Apply group write permissions to all directories above (for collaboration or shared environments)
chmod -R g+w ../*

# Move into the BAM directory to simplify file handling
cd "$bam_dir" || { echo "Error: Cannot change to BAM directory: $bam_dir"; exit 1; }

# Loop over all BAM files in the directory
for bam in *.bam; do
    # Check if the BAM file exists
    if [[ ! -f "$bam" ]]; then
        echo "Warning: No BAM files matching pattern: $bam"
        continue
    fi

    sample=$(basename "$bam" .bam)  # Get sample name by stripping the .bam extension
    echo "Processing sample: $sample"

    # Always reindex the BAM file to ensure accuracy
    echo "Reindexing BAM file: $bam"
    samtools index "$bam" || { echo "Error indexing BAM file: $bam"; continue; }

    # 1. Call variants using FreeBayes (haplotype-aware)
    echo "Calling variants for: $sample using FreeBayes"
    freebayes -f "$ref" \
        --ploidy 2 \
        --min-base-quality 13 \
        --min-mapping-quality 20 \
        --min-coverage 3 \
        --genotype-qualities \
        "$bam" > "${vcf_dir}/${sample}_ace2_raw.tmp" || { echo "Error calling variants with FreeBayes for $sample"; continue; }

    # Check if the raw VCF file was created successfully
    if [[ ! -s "${vcf_dir}/${sample}_ace2_raw.tmp" ]]; then
        echo "Error: Raw VCF file is empty for $sample"
        continue
    fi

    # Rename the temporary file to the raw VCF
    mv "${vcf_dir}/${sample}_ace2_raw.tmp" "${vcf_dir}/${sample}_ace2_raw.vcf"

    # Compress and index the VCF file
    echo "Compressing and indexing ace2_raw VCF: ${sample}_ace2_raw.vcf"
    bgzip -c "${vcf_dir}/${sample}_ace2_raw.vcf" > "${vcf_dir}/${sample}_ace2_raw.vcf.gz" || { echo "Error compressing VCF for $sample"; continue; }
    tabix -p vcf "${vcf_dir}/${sample}_ace2_raw.vcf.gz" || { echo "Error indexing ace2_raw VCF for $sample"; continue; }

    # 4. Filter for depth (DP >= 5) and specify region ### FINAL VERSION CHANGE TO DP >= 10
    echo "Filtering ace2_raw VCF for depth (DP >= 5) and region (ON563187.1:1-712): $sample"  
    bcftools view -r "ON563187.1:1-712" -i 'INFO/DP>=5 || FORMAT/DP>=5' "${vcf_dir}/${sample}_ace2_raw.vcf.gz" -Oz -o "${vcf_dir}/${sample}_filtered_ace2_phased.vcf.gz" || { echo "Error filtering VCF for $sample"; continue; }

    # Check if the filtered VCF file is non-empty
    if [[ ! -s "${vcf_dir}/${sample}_filtered_ace2_phased.vcf.gz" ]]; then
        echo "Warning: Filtered VCF is empty for $sample, skipping consensus generation"
        continue
    fi

    # 5. Index the filtered ace2_phased VCF file
    echo "Indexing filtered ace2_phased VCF: ${sample}_filtered_ace2_phased.vcf.gz"
    tabix -p vcf "${vcf_dir}/${sample}_filtered_ace2_phased.vcf.gz" || { echo "Error indexing filtered VCF for $sample"; continue; }

    # 6. Generate a consensus sequence from the ace2_phased VCF if it contains variants
    if bcftools view -H "${vcf_dir}/${sample}_filtered_ace2_phased.vcf.gz" | head -n 1 | grep -qv "^#"; then
        echo "Generating consensus sequence for: $sample"
        bcftools consensus -f "$ref" "${vcf_dir}/${sample}_filtered_ace2_phased.vcf.gz" > "${con_dir}/${sample}_ace2_phased_consensus.fa" || { echo "Error generating consensus for $sample"; continue; }
        echo "Extracting haplotypes for: $sample"
        vcftools --vcf ${vcf_dir}/${sample}_filtered_ace2_phased.vcf.gz --phased --extract-FORMAT-info GT --out ${hap_dir}/${sample}_filtered_ace2_haplotypes
        awk '{print ">"$1"_"$2"\n"$5"\n>"$1"_"$2"_alt\n"$6}' ${hap_dir}/${sample}_filtered_ace2_haplotypes > "${hap_dir}/${sample}_filtered_ace2_haplotypes.fa" || { echo "Error generating haplotypes for $sample"; continue; }   
    else
        echo "No variants found in filtered VCF for $sample, skipping consensus"
    fi

    # Cleanup intermediate VCF files to save space
    echo "Cleaning up intermediate files for: $sample"
    rm -f "${vcf_dir}/${sample}_ace2_raw.vcf.gz" "${vcf_dir}/${sample}_ace2_raw.vcf.gz.tbi"

    echo "Completed processing for: $sample"
done

echo "All BAM files processed successfully."

# Reapply group write permissions to the full directory (ensures output is accessible)
chmod -R g+w ../*
