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

## Check depth for a single individual, focus on quinx genotypes for Ace2
```
for file in \
B504-UT-M07101-240702_S188_L001.depth.txt \
B053-UT-M07101-240702_S28_L001.depth.txt \
B372-UT-M07101-240702_S105_L001.depth.txt \
B373-UT-M07101-240702_S106_L001.depth.txt \
B374-UT-M07101-240702_S107_L001.depth.txt \
B377-UT-M07101-240702_S108_L001.depth.txt \
B378-UT-M07101-240702_S109_L001.depth.txt \
B379-UT-M07101-240702_S110_L001.depth.txt \
B380-UT-M07101-240702_S111_L001.depth.txt \
B381-UT-M07101-240702_S112_L001.depth.txt \
B393-UT-M70330-240718_S172_L001.depth.txt \
B503-UT-M07101-240702_S187_L001.depth.txt \
B109-UT-M70330-240718_S113_L001.depth.txt \
B263-UT-M70330-240718_S140_L001.depth.txt \
B292-UT-M70330-240718_S154_L001.depth.txt \
B302-UT-M70330-240718_S158_L001.depth.txt \
B392-UT-M70330-240718_S171_L001.depth.txt
do
  gnuplot -persist <<EOF
  set terminal dumb size 120,30
  set title "Read Depth (capped at 20×): ${file}"
  set xlabel "Position (line order)"
  set ylabel "Depth"
  set yrange [0:20]
  plot "${file}" using 0:3 with lines notitle
EOF

  echo "Press enter to continue to the next file..."
  read
done
```
By default, vcftools removes genotypes (entire positions) that don’t pass filters like --minDP. This leads to:

Sites absent from the filtered VCF entirely.

Which bcftools consensus interprets as: “use the reference base.”

So you end up with a consensus that looks complete but isn't accurately masked — especially when there's no data at those sites.
