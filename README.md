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
/uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_raw  

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
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon*
```

# Submit your job

## Github
***Run just ONCE: Clone from Github***  
Before running, you need to make these files on github, and then use git to clone
```
# Just once:
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC              #Example for Emily with "EC"
git clone https://github.com/saarman/cx_amplicon cx_amplicon_scripts
```
***Need to run every time: Pull any changes from Github***  
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_EC/cx_amplicon_scripts/
git pull
```
***Need to run every time to Submit to slurm (change file for each step)***
https://www.chpc.utah.edu/documentation/software/slurm.php#usingslurm
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_scripts
sbatch 1_fastqc.slurm
```

##  Checking the status of your job
https://www.chpc.utah.edu/documentation/software/slurm.php#squeue  
```
squeue --me
```

## Change permissions (again)
```
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_*
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
     - samclip: alignment-aware trimming, remove unwanted primers/adapters from aligned reads (step not in ddRAD)
     - samtools sort: get ready for variant calling
6. bwa flagstat: summarize assembly results (step 6 in ddRAD)
     - sambamba: per base coverage/depth statistics  (optional extra information)
     - samtools idxstats aligned.bam  # summary of how many reads aligned to each genus (optional extra information)
     - samtools view -b -q 30 aligned.bam > high_confidence_alignments.bam # filter out low-confidence mappings (optional extra information)

## Preliminary output
7. vcftools to make vcf file
8. bcftools to call consensus sequence, preliminary analysis only

## Future implementations:
7. freebayes: haplotype-aware variant calling (phasing included, step not in ddRAD)
8. vcftools: filtering snps and haplotypes, or diverge entirely to analyze haplotypes (optional extra information)

# Viewing results files:

## View a BAM file:
```
module load samtools
samtools mpileup -f reference.fasta -Q 20 -q 20 yourfile.bam | less -S
```

## View a VCF file:
```
module load bcftools
bcftools view yourfile.vcf | less -S
```

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
## Norah's
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_bwa  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_fastqc  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_scripts  
 /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_vcf 

 ### Example starting from beginning for Norah:
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS                         # run only once!
git clone https://github.com/saarman/cx_amplicon cx_amplicon_scripts                     # run only once!
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/cx_amplicon_scripts     # run every time
git pull                                                                                 # run every time
sbatch 1_fastqc.slurm                                                                    # run example for step 1
git pull                                                                                 # run every time
sbatch 2_multiqc_summary.slurm                                                           # run example for step 2
```  
# Look at depth stats
Commands:

```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/cx_amplicon_NS/
cd ../cx_amplicon_bwa/depth
awk '$1 == "KY929304.1" && $2 >= 880 && $2 <= 1834' mean_depth_per_position.tsv > mean_depth_KY929304.1_880-1834.tsv

gnuplot -persist <<EOF
set terminal dumb size 120,30
set title "Mean Depth: KY929304.1 (880–1834)"
set xlabel "Position"
set ylabel "Depth"
plot "mean_depth_KY929304.1_880-1834.tsv" using 2:3 with lines title "Mean Depth"
EOF
```
Output:
```

                                                                                                                        
                                                Mean Depth: KY929304.1 (880–1834)                                       
                                                                                                                        
       25000 +------------------------------------------------------------------------------------------------------+   
             |        +         +        +        +         +        +         +        +        +         +        |   
             |                                                                                   Mean Depth ******* |   
             |                                           *** *****                                                  |   
             |                              **       *********   ****                                               |   
       20000 |-+                          *****     ** **           ********                                      +-|   
             |                           ***  **  ***                      **********                               |   
             |                         ***    *** *                          ****  ***                              |   
             |                       ***        ***                                  **                             |   
       15000 |-+                    **          *                                     **                          +-|   
             |                     **                                                  ****                         |   
             |                   ***                                                      **                        |   
             |                 ***                                                         ***                      |   
             |                **                                                             **                     |   
       10000 |-+            ***                                                               *                   +-|   
             |            ***                                                                 ***                   |   
             |          **                                                                      **                  |   
             |         **                                                                        ****               |   
        5000 |-+      **                                                                            **            +-|   
             |        *                                                                              **             |   
             |        *                                                                               **            |   
             |       **                                                                                **           |   
             |      **+         +        +        +         +        +         +        +        +      ** +        |   
           0 +------------------------------------------------------------------------------------------------------+   
            800      900       1000     1100     1200      1300     1400      1500     1600     1700      1800     1900 
                                                            Position                                                    
```



