## RNA_seq differential expression
# Prepare the data
Download the reference genome from Ensembl and use the human GRCh38 version of the genome.

 wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

# Preparing the environement 

 conda create -y --name ngs1 python=3.6

# Data Retrieval
NCBI’s fastq-dump from sra-toolkit was used to download the short reads for NCBI short read archive (SRA). 

 cd ~/workdir/sample_data

 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 -N 10000 -X 4010000    --clip SRR1039520


# Data Retrieval(past troubleshooting)

Direct download from NCBI

  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz

#Using SRA-toolkit 

  Prefetch SRR103950
  fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039508
  prefetch SRR1039509
  fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039509
 
 
 The data was quite large around 3 gigabyte each and consuming more than 3 hours download time let alone editing it 

Trying to get subset from the net 
Wget airway package 

 # Setup enviornemnt

conda activate ngs1

conda install -c bioconda fastqc 

conda install -c bioconda multiqc

conda install sra-toolkit\

#conda install samtools

conda install -c bioconda -y hisat2

#conda install kallisto

# install r and dependicies
conda install r
conda install -y bioconductor-deseq r-gplots

 

Got the bam files 



# Analyzing control samples
# Alignment: Hisat2

# Step 1 (Indexing)

INDEX=./gencode.v33.transcripts

REF= ./gencode.v33.transcripts.fa


hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv gencode.v33.transcripts.fa gencode.v33.transcripts
 REF_ERCC=./ref/ERCC92.fa

hisat2-build $REF $INDEX


# Step 2 (Alignment)

RUNLOG=runlog.txt
READS_DIR=~/Downloads/fastqq/fastq/fastg/
mkdir bam
# align both reads for treated and untreated conditions

for SAMPLE in UNT;
do
    for REPLICATE in 12 16 20;
    do
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq.gz
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq.gz
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done

Check the Alignment summary.
4000001 reads; of these:
  4000001 (100.00%) were paired; of these:
    656023 (16.40%) aligned concordantly 0 times
    894885 (22.37%) aligned concordantly exactly 1 time
    2449093 (61.23%) aligned concordantly >1 times
    ----
    656023 pairs aligned concordantly 0 times; of these:
      6338 (0.97%) aligned discordantly 1 time
    ----
    649685 pairs aligned 0 times concordantly or discordantly; of these:
      1299370 mates make up the pairs; of these:
        1107153 (85.21%) aligned 0 times
        61740 (4.75%) aligned exactly 1 time
        130477 (10.04%) aligned >1 times
86.16% overall alignment rate

for SAMPLE in TTT;
do
    for REPLICATE in 13 17 21;
    do
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq.gz
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq.gz
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done

Check the Alignment summary.
4000001 reads; of these:
  4000001 (100.00%) were paired; of these:
    644104 (16.10%) aligned concordantly 0 times
    925147 (23.13%) aligned concordantly exactly 1 time
    2430750 (60.77%) aligned concordantly >1 times
    ----
    644104 pairs aligned concordantly 0 times; of these:
      6320 (0.98%) aligned discordantly 1 time
    ----
    637784 pairs aligned 0 times concordantly or discordantly; of these:
      1275568 mates make up the pairs; of these:
        1061485 (83.22%) aligned 0 times
        66172 (5.19%) aligned exactly 1 time
        147911 (11.60%) aligned >1 times
86.73% overall alignment rate




# Step 3 (Quantification)

GTF=~/Downloads/fastqq/fastq/fastg/gencode.v33.transcripts.annotation.gtf 


# Generate the counts.
featureCounts -a $GTF -g gene_name -o counts.txt  bam/UNT*.bam  bam/TTThe chromosome name of "ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|" contains 125 characters, longer than the upper limit of 99
featureCounts has to stop runningT*.bam

# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt
#Head
https://github.com/ngs-project/RNA_seq/blob/master/simple_counts.txt


The chromosome name of "ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|" contains 125 characters, longer than the upper limit of 99
featureCounts has to stop running
||                                                                            ||
||    format error found in this file!                                        ||

FATAL Error: The program has to terminate and no counting file is generated.



# Analyze the counts with DESeq1.
cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
https://github.com/ngs-project/RNA_seq/blob/master/results_deseq1.tsv

DeSEQ1 Output header description

    id: Gene or transcript name that the differential expression is computed for,
    baseMean: The average normalized value across all samples,
    baseMeanA, baseMeanB: The average normalized gene expression for each condition,
    foldChange: The ratio baseMeanB/baseMeanA ,
    log2FoldChange: log2 transform of foldChange . When we apply a 2-based logarithm the values become symmetrical around 0. A log2 fold change of 1 means a doubling of the expression level, a log2 fold change of -1 shows show a halving of the expression level.
    pval: The probability that this effect is observed by chance,
    padj: The adjusted probability that this effect is observed by chance.

View only rows with pval < 0.05

cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf














