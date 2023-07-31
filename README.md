# RNAseq-pancreatic-cells

# Overview
This project is a replication of an RNA seq analysis of lncRNA’s in human pancreatic islets and β-cells, with the goal of uncovering interesting insights about the expression of those cells that correspond with those of the original study. The analysis was performed in Linux Bash utilizing famous cli tools used in RNA seq analyses. An extensive review of the operations performed and the corresponding results is presented in the Report file. 

# Datasets
The original study from which the data were derived was the following: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3475176/ . The data were downloaded from "ebi.ac.uk" and were were FASTQ files from a transcriptome analysis of lncRNA’s in human pancreatic islets and β-cells. There were four files in total, two for each sample. They were paired end data with each part of the file being the left or right read.  The original study defined 1128 islet lncRNA genes, showing them to be an integral component of the dynamic β-cell specific differentiation program.

# Methods
A typical RNA seq analysis scheme was followed for most of the analysis. The operations performed (in order of performance) were:
* Quality Control and Preprocessing
    * Extraction of basic sequence statistics
    * Sequence quality per base
    * Sequence quality per file
    * Quality scores per sequence
    * Sequence content per base
    * GC content per sequence
    * N content per base
    * Sequence length distribution
    * Sequence duplication levels
    * Overrepresented sequences
    * Adapter content

* Identification of unknown adapters with the Minion tool, developed by the Genomic's Lab in Fleming Research Institute (Athens, Greece)
* Removal of unwanted and overrepresented sequences with Cutadapt
* Spliced alignment (hg38 used as reference genome)
* Differential expression analysis with cutdiff (from Cufflinks package)
* Visualization of the top 5 differentially expressed genes with the Integrative Genomics Viewer (IGV) tool
* GO Term enrichment analysis for the significant genes with g::Profiler

# Results
The five highest differentially expressed genes were LGI1, MIR210HG, SLC37A4, SERPINA3 and MIR492. The gene LGI1 is involved in the development and morphogenesis of neuronal projections, MIR210HG plays a role in the cellular response to hypoxia and angiogenesis, SLC37A4 is involved in glucose and carbohydrate transport, SERPINA3 is involved in the acute-phase and inflammatory response, and MIR492 is involved in cell proliferation and transcriptional regulation. 
The seven significantly GO terms that were identified in the GO terms enrichment analysis were homeostatic process, regulation of biological quality, response to stimulus, regulation of signaling, regulation of type IIa hypersensitivity, regulation of cell communication, apical plasma membrane and extracellular region. The significant GO terms suggest that the biological processes that the significantly expressed genes are involved with are regulated and maintained by homeostatic processes, response to stimuli, regulation of signaling and communication, and the extracellular environment.
Based on these results of the transcriptome analysis of lncRNAs in human pancreatic islets and β-cells, it can be infered that these cells play important roles in regulating glucose and carbohydrate transport, acute-phase and inflammatory response, and cellular response to hypoxia and angiogenesis. The significant GO
terms identified provide additional insights into the potential regulatory pathways involved in the maintenance and function of these cells, which are consistent with their known functions.








