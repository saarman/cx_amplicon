# How Reads Are Aligned to Multiple References
If multiple references are provided in the FASTA file, each read is aligned to all of them.
BWA MEM assigns **each read to the reference with the highest alignment score

# Handling Multiple Alignments
Primary Alignment (Default Behavior):

By default, only the best-scoring alignment is reported.
This means each read is assigned to one location in the reference genome.
Secondary Alignments (-a option):

If a read has multiple high-scoring alignments, you can ask bwa mem to report secondary alignments using:
```
bwa mem -a reference.fasta reads.fastq > alignments.sam
```

# Secondary alignments appear in the SAM file with the 0x100 flag.
Supplementary Alignments (-M option):

For split-read alignments (e.g., structural variants or chimeric reads), the -M flag marks supplementary alignments.
Unmapped Reads:

Reads that do not align well anywhere are marked as unmapped (* in the reference column of the SAM file).


### Creating consensus for each sample and adding name of sample to header:
https://chatgpt.com/share/67e4403c-ab24-800e-b5e6-8e0653b060e6

```
for bam in *.bam; do
  sample=$(basename "$bam" .bam)

  # Ensure BAM is indexed
  if [ ! -f "${bam}.bai" ] && [ ! -f "${sample}.bai" ]; then
    samtools index "$bam"
  fi

  echo "Processing $sample..."

  samtools mpileup -uf reference.fasta "$bam" | \
  bcftools call -c | \
  bcftools filter -s LowQual -e '%QUAL<20 || DP<10' | \
  bcftools view -R scaffolds.txt | \
  bcftools consensus -f reference.fasta | \
  awk -v sample="$sample" '/^>/{print ">"sample"|"substr($0,2)} !/^>/' > "${sample}_consensus.fa"
done
```

Still some tweaks...testing that input and output exists...
```

bash

# Load required software modules
module load bwa/2020_03_19
module load samtools/1.16
module load bcftools/1.16

# Define input/output paths and filenames
bam_dir="./../cx_amplicon_bwa"                    # Directory containing input BAM files
vcf_dir="./../cx_amplicon_vcf"                    # Directory to store per-sample VCF files
out_dir="./../cx_amplicon_consensus"              # Directory to store consensus FASTA files
ref="../cx_amplicon_bwa/ref/Rep_Genera_Mito.fasta" # Reference FASTA file
region_file="ref_list_ace2.txt"                   # Text file with one region: e.g., ON563187.1:1-712
gene="ace2"                                       # Target gene name for file labeling

# Create output directories if they don't exist
mkdir -p "$out_dir"
mkdir -p "$vcf_dir"

# Ensure group write permissions for relevant directories
chmod -R g+w "$out_dir" "$vcf_dir"

samtools faidx "$ref" ON563187.1:1-712

```
## Ideas for improvement
Separate pipelines for COI from mosquito genes, from step 3 on.
COI: use a consensus for vertebrates as the single reference, and use extremely lenient settings
Mosquito: allow for gap opening, especially for cqm1 
