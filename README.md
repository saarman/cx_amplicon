# cx_amplicon
Culex species identification, cqm1 analysis, and blood meal identification with bwa mem assembly of illumina reads from amplicon resequencing (COI, Ace2, cqm1)

## A place to keep track of amplicon reference assembly bioinformatics steps

# Logging onto CHPC with Terminal on a mac
1. Open Terminal
2. ssh u6036559@notchpeak.chpc.utah.edu        #replace with your username
3. salloc --time=72:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np

# Logging onto CHPC with Remote Desktop interactive job
1. Go to "my interactive sessions" page on CHPC [CHPC](https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sessions)
2. In the lefthand menu, choose "Interactive Desktop"
3. Fill in dropdown menu with this info:
   - Cluster: notchpeak
   - Account, partition, qos: saarman-np:saarman-shared-np:saarman-np
   - Number of cores (per node): 1
   - Number of hours: 72             # up to 336 allowed
   - Advanced options, mem per job: 10       # up to 100 gb

## Upload files to CHPC group1 storage
SFTP file transfer protocol with Cyberduck worked easily.  
Select SFTP at the top.  
Server: notchpeak  
User: u6036559  
Port: 21  
Drag and drop from local disk.  
/uufs/chpc.utah.edu/common/home/saarman-group1/cx_ddRAD_raw  

## To unzip a .tar file (Norah completed 11/8/2024):
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_ddRAD_raw
tar -xvf Culex_ddRAD.tar 
mv ./Culex2/Culex2.fastq.gz cx_ddRAD_plate2.fastq.gz
```

## To unzip a .tar.gz file:
```
tar –xvzf cx_ddRAD_plate2.tar.gz
mv ./data/Saarman/ddRAD/ddRAD.fastq.gz ./cx_ddRAD_plate2.fastq.gz
```

## Change permissions
```
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_ddRAD*
```

# Submit your job

## Github
***Run just ONCE: Clone from Github***  
Before running, you need to make these files on github, and then use git to clone
```
# Just once:
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC              #Example for Emily with "EC"
git clone https://github.com/saarman/cx_amplicon_scripts cx_amplicon_scripts
```
***Need to run every time: Pull any changes from Github***  
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_scripts/
git pull
```
***Need to run every time: Submit to slurm***
https://www.chpc.utah.edu/documentation/software/slurm.php#usingslurm
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_ddRAD_scripts/
sbatch 1a_process_radtags.slurm
```

##  Checking the status of your job
https://www.chpc.utah.edu/documentation/software/slurm.php#squeue  
```
squeue --me
```

## Change permissions (again)
```
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_ddRAD*
```

# Reference files 

Ref files are saved as _REF.fasta
The final file we plan to use for cqm1 and ace2 workflow is:  
culex_cqm1_ace2_REF.fasta

# Pipeline steps 
 
1. fastqc: Read quality (step 2 in ddRAD)
2. multiqc: Summarize quality results (step 3 in ddRAD)
4. bwa index: build bwa index for reference fasta files (step 4b in ddRAD)
5. bwa mem: assemble to reference file with .pl and .slurm commands (step 5 in ddRAD) 
     - samclip: alignment-aware trimming, specifically designed to remove unwanted sequences like primers or adapters from aligned reads (step not in ddRAD)
     - samtools sort: get ready for variant calling
6. bwa flagstat: summarize assembly results (step 6 in ddRAD)
     - sambamba: per base coverage/depth statistics  (optional extra information)
     - samtools idxstats aligned.bam  # summary of how many reads aligned to each genus (optional extra information)
     - samtools view -b -q 30 aligned.bam > high_confidence_alignments.bam # filter out low-confidence mappings (optional extra information)
7. freebayes: haplotype-aware variant calling (phasing included, step not in ddRAD)
8. vcftools: filtering snps and haplotypes, or diverge entirely to analyze haplotypes (optional extra information)


# Directory structure:

## Raw files from uphl .fastq.gz
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_raw                     

 ## Note that references will be on github 
 So will be pulled down into the cx_amplicon_scripts folders

 ## Emily's 
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC   
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_bwa  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_fastqc  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_scripts  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_vcf  
 ## Eric's
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EJ   
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EJ/cx_amplicon_bwa  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EJ/cx_amplicon_fastqc  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EJ/cx_amplicon_scripts  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EJ/cx_amplicon_vcf  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_vcf  
 ## Matt's 
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_MY  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_MY/cx_amplicon_bwa  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_MY/cx_amplicon_fastqc  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_MY/cx_amplicon_scripts  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_MY/cx_amplicon_vcf  
 ## Katie's  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_KG  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_KG/cx_amplicon_bwa  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_KG/cx_amplicon_fastqc  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_KG/cx_amplicon_scripts  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_KG/cx_amplicon_vcf  
 ## Tyler's
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_TS  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_TS/cx_amplicon_bwa  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_TS/cx_amplicon_fastqc  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_TS/cx_amplicon_scripts  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_TS/cx_amplicon_vcf  
