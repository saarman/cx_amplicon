#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;

my $max = 20;  # Set the maximum number of parallel processes to 1 for testing, target is 20
my $pm = Parallel::ForkManager->new($max);  # Create a new Parallel::ForkManager object with the specified maximum

# Path to the reference file
my $ref = "/uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_scripts/Rep_Genera_Mito.fasta";

# Path to the primers file
my $primers = "/uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_scripts/primers.fasta";

# Output directory
my $output_dir = "/uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_bwa";

# Path to samtools
my $samtools = "/uufs/chpc.utah.edu/sys/installdir/samtools/1.16/bin/samtools"; # module load samtools/1.16; which samtools

# Path to bwa-mem2 binary, find this with: 
my $bwa = "/uufs/chpc.utah.edu/sys/installdir/bwa/2020_03_19/bin/bwa"; # module load bwa/2020_03_19; which bwa

FILES:
foreach my $fq1 (@ARGV) {  # Iterate over each file passed as an argument
    $pm->start and next FILES;  # Fork a new process and move to the next file if in the parent process

    # Extract the identifier from the filename
    $fq1 =~ m/([A-Za-z_\-0-9]+)\.fq\.gz$/ or die "failed match for file $fq1\n";
    my $ind = $1;  # Store the identifier in $ind

    # Run the BWA-MEM2 alignment and process with samtools, could add -K 1000000 -c 1000 to reduce mem?
    my $cmd = "$bwa mem -M -t 1 $ref $fq1 | samclip --ref $primers --max 50 | $samtools view -b | $samtools sort --threads 1 > ${output_dir}/${ind}.bam";
    system($cmd) == 0 or die "system $cmd failed: $?";

    print "Alignment completed for $ind\n";

    $pm->finish;  # End the child process
}

$pm->wait_all_children;  # Wait for all child processes to finish





# bwa-mem2 mem #options
#-----------------------------
#Executing in AVX512 mode!!
#-----------------------------
#Usage: bwa2 mem [options] <idxbase> <in1.fq> [in2.fq]
#Options:
#  Algorithm options:
#    -o STR        Output SAM file name
#    -t INT        number of threads [1]
#    -k INT        minimum seed length [19]
#    -w INT        band width for banded alignment [100]
#    -d INT        off-diagonal X-dropoff [100]
#    -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
#    -y INT        seed occurrence for the 3rd round seeding [20]
#    -c INT        skip seeds with more than INT occurrences [500]
#    -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
#    -W INT        discard a chain if seeded bases shorter than INT [0]
#    -m INT        perform at most INT rounds of mate rescues for each read [50]
#    -S            skip mate rescue
#    -o            output file name missing
#    -P            skip pairing; mate rescue performed unless -S also in use
#Scoring options:
#   -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
#   -B INT        penalty for a mismatch [4]
#   -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
#   -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
#   -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
#   -U INT        penalty for an unpaired read pair [17]
#Input/output options:
#   -p            smart pairing (ignoring in2.fq)
#   -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
#   -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
#   -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)
#   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
#   -T INT        minimum score to output [30]
#   -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
#   -a            output all alignments for SE or unpaired PE
#   -C            append FASTA/FASTQ comment to SAM output
#   -V            output the reference FASTA header in the XR tag
#   -Y            use soft clipping for supplementary alignments
#   -M            mark shorter split hits as secondary
#   -I FLOAT[,FLOAT[,INT[,INT]]]
#                 specify the mean, standard deviation (10% of the mean if absent), max
#                 (4 sigma from the mean if absent) and min of the insert size distribution.
#                 FR orientation only. [inferred]
#Note: Please read the man page for detailed description of the command line and options.
#No. of OMP threads: 0
#Processor is runnig @2095.190746 MHz
#Runtime profile:
#
#	 Time taken for main_mem function: 0.00 Sec
#
#	IO times (sec) :
#	Reading IO time (reads) avg: 0.00, (0.00, 0.00)
#	Writing IO time (SAM) avg: 0.00, (0.00, 0.00)
#	Reading IO time (Reference Genome) avg: 0.00, (0.00, 0.00)
#	Index read time avg: 0.00, (0.00, 0.00)
#
#	Overall time (sec) (Excluding Index reading time):
#	PROCESS() (Total compute time + (read + SAM) IO time) : 0.00
#	MEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: 0.00, (0.00, 0.00)
#
#	 SAM Processing time (sec):
#	--WORKER_SAM avg: 0.00, (0.00, 0.00)
#
#	Kernels' compute time (sec):
#	Total kernel (smem+sal+bsw) time avg: 0.00, (0.00, 0.00)
#		SMEM compute avg: -nan, (0.00, 477283.51)
#		SAL compute avg: -nan, (0.00, 477283.51)
#		BSW time, avg: -nan, (0.00, 477283.51)
#
#	Total allocs: 0 = 0 out total requests: 0, Rate: -nan
#
#Important parameter settings: 
#	BATCH_SIZE: 512
#	MAX_SEQ_LEN_REF: 256
#	MAX_SEQ_LEN_QER: 128
#	MAX_SEQ_LEN8: 128
#	SEEDS_PER_READ: 500
#	SIMD_WIDTH8 X: 64
#	SIMD_WIDTH16 X: 32
#	AVG_SEEDS_PER_READ: 64#
