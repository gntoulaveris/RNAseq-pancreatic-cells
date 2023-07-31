# This is a transcriptome analysis of lncRNA’s in human pancreatic islets and β-cells.

# Data acquirement
# the same command for all files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR173/ERR173261/ERR173261_2.fastq.gz
# unzip them all
for file in *.fastq.gz; do gunzip "$file"; done

# Quality Control and Preprocessing
cd FastQC

unzip fastqc_v0.11.9.zip

# grant execute permissions to the file
chmod +x fastqc

# run file
./fastqc


# Identify unknown parameters in the samples - Minion
cd ~/tools
wget https://genomics-lab.fleming.gr/fleming/uoa/vm/minion
chmod u+x minion
~/tools/minion --help

~/tools/minion search-adapter -i ERR173261_1.fastq
~/tools/minion search-adapter -i ERR173261_2.fastq

~/tools/minion search-adapter -i ERR173280_1.fastq
~/tools/minion search-adapter -i ERR173280_2.fastq


# Remove unwanted overrepresented seqs - Cutadapt
# install pip
tar xvfz pip-23.0.tar.gz
cd pip-23.0
python setup.py install --user
# make sure the installation was successful
/home/grad/.local/bin/pip3.8

# pip was added to bash PATH so that it could be called from any directory
nano ~/.bashrc # open bash configuration file
export PATH=$PATH:/home/grad/.local/bin # written at the end of the file

# reload terminal to set the change in effect
source ~/.bashrc

# setuptools were requested to install cutadapt
pip install setuptools_scm

pip install cutadapt

cutadapt -u 30 -o ERR173261_2_trimmed.fastq ERR173261_2.fastq


# Spliced alignment
tar xvfz genome_assemblies_genome_gtf.tar
cd ncbi-genomes-2023-02-06
gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
mv GCF_000001405.40_GRCh38.p14_genomic.gtf hg38.gtf
mv GCF_000001405.40_GRCh38.p14_genomic.gtf ~/  #move to home

tar -xvf genome_assemblies_genome_fasta.tar
cd ncbi-genomes-2023-02-08
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna hg38.fasta
mv hg38.fasta ~/  #move to home

cd ~/


# a bowtie index file of hg38 must be created first
bowtie2-build hg38.fa hg38

# alignment with tophat
tophat2 -p 8 --no-coverage-search -o ERR173261_aligned hg38 ERR173261_1.fastq.gz ERR173261_2.fastq.gz

tophat2 -p 8 --no-coverage-search -o ERR173280_aligned hg38 ERR173280_1.fastq.gz ERR173280_1.fastq.gz


# Differential expression analysis
```bash
# a bowtie index file of hg38 must be created first
bowtie2-build hg38.fa hg38

# alignment with tophat
tophat2 -p 8 --no-coverage-search -o ERR173261_aligned hg38 ERR173261_1.fastq.gz ERR173261_2.fastq.gz

tophat2 -p 8 --no-coverage-search -o ERR173280_aligned hg38 ERR173280_1.fastq.gz ERR173280_1.fastq.gz

grep "yes" "gene_exp.diff" | grep -v "^test_id" > significant_diff_expressed_genes.txt

test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant

grep "yes" "gene_exp.diff" | grep -v "^test_id" | cut -f2,3


# IGV Visualization of top 5 differentially expressed genes
grep "yes" gene_exp.diff | grep -v "^test_id" | sort -k 12n,12 -k 11nr,11 | head -n 15 | cut -f 2,3

# create the index bai files
samtools index ERR173261_accepted_hits.bam

samtools index ERR173280_accepted_hits.bam

# check if BAM files are properly aligned
samtools view -H ERR173261_accepted_hits.bam | grep -w "SO:coordinate"
@HD	VN:1.0	SO:coordinate

samtools view -H ERR173280_accepted_hits.bam | grep -w "SO:coordinate"
@HD	VN:1.0	SO:coordinate

samtools view -h ERR173261_accepted_hits.bam | head -30

# Open IGV

# GO term enrichment of significant genes
cut -f1,2 significant_diff_expressed_genes.txt > filtered_genes.txt




















