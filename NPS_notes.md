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
