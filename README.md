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

  Prefetch SRR1039508

Then fast_dump
  fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039508
  prefetch SRR1039509
  fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039509
 
 
 The data was quite large around 3 gigabyte each and consuming more than 3 hours download time let alone editing it 

Trying to get subset from the net 
Wget airway package 



 # Installation 

Conda activate ngs1

conda activate ngs1

conda install -c bioconda fastqc 

conda install -c bioconda multiqc

Conda install sra-toolkit\

Conda install samtools

Conda install -c bioconda -y hisat2

Conda install kallisto

 

Got the bam files 
