## RNA_seq differential expression
# Prepare the data
Download the reference genome from Ensembl and use the human GRCh38 version of the genome.

 wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

# Preparing the environement 

 conda create -y --name ngs1 python=3.6

# Data Retrieval

 cd ~/workdir/sample_data

 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 -N 10000 -X 4010000    --clip SRR1039520

NCBI’s fastq-dump from sra-toolkit was used to download the short reads for NCBI short read archive (SRA). 

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

Conda activate ngs1

conda install -c bioconda fastqc 

conda install -c bioconda multiqc

Conda install sra-toolkit\

#Conda install samtools

Conda install -c bioconda -y hisat2

#Conda install kallisto

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

# Step 3 (Quantification)

GTF=~/Downloads/fastqq/fastq/fastg/

# Generate the counts.
featureCounts -a $GTF -g gene_name -o counts.txt  bam/TTT*.bam  bam/UNT*.bam

# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt
#Head














