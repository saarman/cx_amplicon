#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;

my $max = 20;  # Set the maximum number of parallel processes
my $pm = Parallel::ForkManager->new($max);  # Create a new Parallel::ForkManager object

# Path to the reference file
my $ref = "../cx_amplicon_bwa/ref/ace2_cqm1_coi.fasta";

# Path to the primers file
my $primers = "../cx_amplicon_bwa/ref/primers.fasta";

# Output directory
my $output_dir = "../cx_amplicon_bwa";

# Path to samtools
my $samtools = "/uufs/chpc.utah.edu/sys/installdir/samtools/1.16/bin/samtools";

# Path to bwa-mem2 binary
my $bwa = "/uufs/chpc.utah.edu/sys/installdir/bwa/2020_03_19/bin/bwa";

# Path to samclip 
my $samclip = "/uufs/chpc.utah.edu/common/home/saarman-group1/samclip/samclip";

FILES:
foreach my $fq1 (@ARGV) {  # Iterate over each file passed as an argument

    # Extract the identifier from the filename
    if ($fq1 =~ m/([A-Za-z_\-0-9]+)_R1_001\.fastq\.gz$/) {
        my $ind = $1;  # Store the identifier in $ind
        my $fq2 = "../../cx_amplicon_raw/${ind}_R2_001.fastq.gz";  # Construct the R2 filename

        # Check if the paired-end file exists
        unless (-e $fq2) {
            warn "Paired file not found: $fq2\n";
            next;
        }

        # Fork a process for parallel execution
        $pm->start and next;

        # my $cmd = "$bwa mem -M -t 4 $ref $fq1 $fq2 | $samclip --ref $primers --max 50 | $samtools view -b | $samtools sort --threads 4 > ${output_dir}/${ind}.bam";

        my $cmd = "$bwa mem -M -t 4 -B 2 -O 4,4 -E 2,2 -k 15 -T 20 $ref $fq1 $fq2 | $samtools view -b | $samtools sort --threads 4 > ${output_dir}/${ind}.bam";
        
        system($cmd) == 0 or die "system $cmd failed: $?";   

        print "Alignment completed for $ind\n";

        $pm->finish;  # End the child process
    } else {
        die "Failed match for file $fq1\n";
    }
}

$pm->wait_all_children;  # Wait for all child processes to finish

###############################################
# Explanation of parameter choices April 1, 2025
###############################################
# -B 2: Lower mismatch penalty, allowing for more mismatches.
# -O 4,4: Lower gap open penalties for both insertions and deletions.
# -E 2,2: Reduced gap extension penalties, making long gaps less penalizing.
# -k 15: Allows shorter seeds, increasing alignment sensitivity.
# -T 20: Lowers the threshold for reporting alignments, capturing more divergent reads.
