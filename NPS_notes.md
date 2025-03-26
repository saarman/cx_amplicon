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
```
done
