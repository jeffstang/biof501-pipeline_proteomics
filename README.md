# BIOF 501 Term Project: RNASeq workflow for differential expression and pathway enrichment 

## Table of Contents:
- [Background and Rationale](#background-and-rationale)
    - [Objective and Aims](#objective-and-aims)
    - [Package Dependencies](#package-dependencies)
    - [Workflow Overview](#workflow-overview)
- [Usage](#usage)
    - [Installation](#installation)
    - [Project Structure and Overview](#project-structure-and-overview)
    - [Quickstart](#quickstart)
    - [Step 1: Preprocess FASTQ files](#step-1-preprocess-fastq-files)
    - [Step 2: Quantify reads using Salmon](#step-2-quantify-reads)
    - [Step 3: Perform differential expression analysis](#step-3-perform-differential-expression-analysis)
    - [Step 4: Perform pathway enrichment analysis](#step-4-perform-pathway-enrichment-analysis)
- [References](#references)

## Background and Rationale
Cellular communication occurs through ligands binding to surface receptors, influencing signaling pathways that trigger phenotypic responses which reflect on cell states and/or disease progression [[1](#references), [2](#references)]. Research on ligand-receptor interactions (LRIs) has gained momentum with single-cell (sc-) omics, which provide high-resolution insights into individual cell states [[3](#references)]. However, the volume of bulk transcriptomic datasets, complemented by extensive clinical data available, compared to sc-omics provide a valuable yet under-utilized resource [[4](#references)]. Thus tooling using bulk transcriptomics to characterize cellular networks can complement single-cell studies [[5](#references)]. 

As omics technologies advance, complexity of data generated increase [[6](#references)]. In combination with complexity, rapidly-evolving terminologies, data quality and reliability, and lack of supporting documentation further create variability and barriers for reproducibility when leveraging publicly available archival data sets [[7](#references)]. 

### Objective and Aims:
A popular method in transcriptomics is the gene set enrichment analysis (GSEA) approach which enables identification of activation or repression of gene sets that share common biology and catalogued into molecular pathways [[8](#references)]. With more LRI databases emerging [[9](#references), [10](#references), [11](#references), [12](#references)], curation and validation of these potential annotations become important. We hypothesize that by repurposing GSEA for LRI annotations, we can validate discovery-based LRI methodologies.
  
**Aim 1 -** To provide an end-to-end workflow solution for analyzing changes in gene expression starting from raw **fastq** files.

**Aim 2 -** To explore gene set enrichment analysis in the context of LRIs described by pre-existing database annotations. A **key assumption** of note is that receptor gene expression changes can be indicative of ligand binding and sufficiently describes the ligand activity (in the context of fgsea, the "pathway").  

### Workflow Overview
The workflow includes the following steps:
1. Retrieve reference files from ENCODE
2. Preprocess raw fastq files - provides quality metrics from `fastqc` [[13](#references)] pre- and post-trimming and trims reads using `trimmomatic` [[14](#references)]
3. Count transcripts and generate count matrix using `salmon` [[15](#references)]
4. Preprocess and analyze differentially expressed genes (DEGs) using `limma-voom` [[16](#references)] 
5. We leverage `fgsea` [[17](#references)] and the example database to provide insight on ligand activity which can reflect on signalling pathway activity

Below is a workflow diagram:

```mermaid
graph TD
    %% Main Workflow
    subgraph MAIN_WORKFLOW ["Main Workflow"]
        direction TB
        RAW_FASTQ[Input: RAW FASTQ files] --> TRIM[TRIM_READS: Trim reads using Trimmomatic]
        TRIM --> SALMON_QUANT[SALMON_QUANT: Quantify transcripts using Salmon]
        SALMON_QUANT --> TXIMPORT[TXIMPORT_PROCESS: Merge all salmon results at the gene-level to produce a counts matrix]
        TXIMPORT --> DEA[LIMMA_VOOM_DEA: Run differential expression analysis]
        DEA --> VOLCANO[ENHANCED_VOLCANO_PLOT: Generate a volcano plot to show DEGs]
        DEA --> PATHWAY[PATHWAY_ENRICHMENT: Perform pathway enrichment analysis and get a rank list of ligands represented by LRI-annotated differentially expressed genes]
    end

    %% FASTQC Workflow
    subgraph FASTQC_WORKFLOW ["FASTQC Workflow"]
        direction TB
        RAW_FASTQ --> FASTQC_RAW[FASTQC_RAW: Generate FASTQC for raw reads] 
        FASTQC_RAW --> TRIM
        TRIM --> FASTQC_TRIMMED[FASTQC_TRIMMED: Generate FASTQC for trimmed reads]
    end

    %% SALMON Workflow
    subgraph SALMON_WORKFLOW ["Auxiliary File Requirements"]
        direction TB
        DOWNLOAD_FASTA[DOWNLOAD_FASTA: Download FASTA] --> SALMON_INDEX[SALMON_INDEX: Index FASTA file]
        DOWNLOAD_GTF[DOWNLOAD_GTF: Download GTF] --> CREATE_TX2GENE[CREATE_TX2GENE: Create tx2gene mapping]
        SALMON_INDEX --> SALMON_QUANT
        DOWNLOAD_GTF --> SALMON_QUANT
        CREATE_TX2GENE --> TXIMPORT
    end

    %% Styling
    style MAIN_WORKFLOW fill:#E3F2FD,stroke:#42A5F5,stroke-width:2px
    style FASTQC_WORKFLOW fill:#FFEBEE,stroke:#E57373,stroke-width:2px
    style SALMON_WORKFLOW fill:#FFF3E0,stroke:#FFA726,stroke-width:2px
    
    %% Input and Styling
    style RAW_FASTQ fill:#C8E6C9,stroke:#388E3C,stroke-width:2px
    style FASTQC_RAW fill:#BBDEFB,stroke:#1976D2,stroke-width:2px
    style FASTQC_TRIMMED fill:#BBDEFB,stroke:#1976D2,stroke-width:2px
    style TRIM fill:#FFCCBC,stroke:#E64A19,stroke-width:2px
    style SALMON_INDEX fill:#FFD54F,stroke:#F57F17,stroke-width:2px
    style SALMON_QUANT fill:#FFD54F,stroke:#F57F17,stroke-width:2px
    style DOWNLOAD_FASTA fill:#FFAB91,stroke:#D84315,stroke-width:2px
    style DOWNLOAD_GTF fill:#FFAB91,stroke:#D84315,stroke-width:2px
    style CREATE_TX2GENE fill:#FFAB91,stroke:#D84315,stroke-width:2px
    style TXIMPORT fill:#90CAF9,stroke:#1E88E5,stroke-width:2px
    style DEA fill:#90CAF9,stroke:#1E88E5,stroke-width:2px
    style VOLCANO fill:#CE93D8,stroke:#8E24AA,stroke-width:2px
    style PATHWAY fill:#CE93D8,stroke:#8E24AA,stroke-width:2px
```

### Package Dependencies
All tools used in the workflow have specified docker containers that are **automatically pulled when launching the workflow**. Please see below for tools, version tags, and the docker image used. 
#### List of Docker Images used (see[ nextflow.config](https://github.com/jeffstang/biof501-term_project/blob/main/nextflow.config)):
```
edgeR 4.0.2
quay.io/biocontainers/bioconductor-edger:4.0.2--r43hf17093f_0

EnhancedVolcano 1.20.0 
quay.io/biocontainers/bioconductor-enhancedvolcano:1.20.0--r43hdfd78af_0

FastQC 0.12.1
staphb/fastqc:0.12.1

fgsea 1.28.0
quay.io/biocontainers/bioconductor-fgsea:1.28.0--r43hf17093f_1 

salmon 1.10.3
combinelab/salmon:1.10.3

trimmomatic 0.36
quay.io/biocontainers/trimmomatic:0.36--4                 

tximport 1.26.0
quay.io/biocontainers/bioconductor-tximport:1.26.0--r42hdfd78af_0
```

## Usage
### Installation
The following versions of software were installed following the [Nextflow Docs](https://www.nextflow.io/docs/latest/install.html)
```
bash-5.2.21
docker-26.1.1
git-2.43.0
nextflow-24.10.0
open-jdk-17.0.10
```
If initial setup requirements are met, please clone this GitHub Repository:
```
git clone https://github.com/jeffstang/biof501-term_project.git
```

### Project Structure and Overview
Upon cloning, the repository should come with a set of base files.
<details>
    <summary> Click here to see a tree overview of the overall directory structure </summary>

```bash
├── README.md
├── bin
│   ├── ligand_enrichment_analysis.R
│   ├── limma_voom.R
│   ├── plot_enhancedVolcano.R
│   └── tximport.R
├── data
│   ├── raw
│   │   ├── SRR24360639.sub_1.fastq.gz
│   │   ├── SRR24360639.sub_2.fastq.gz
│   │   ├── SRR24360643.sub_1.fastq.gz
│   │   ├── SRR24360643.sub_2.fastq.gz
│   │   ├── SRR24360647.sub_1.fastq.gz
│   │   ├── SRR24360647.sub_2.fastq.gz
│   │   ├── SRR24360653.sub_1.fastq.gz
│   │   └── SRR24360653.sub_2.fastq.gz
│   └── reference
│       ├── CellChatDB_preprocess.R
│       ├── cellchatv2_mouseLRI.rda
│       └── metadata.csv
├── main.nf
├── modules
│   ├── convert_tx2gene
│   │   └── main.nf
│   ├── diff_exp_analysis
│   │   └── main.nf
│   ├── download_ref
│   │   └── main.nf
│   ├── enhanced_volcano
│   │   └── main.nf
│   ├── fastqc
│   │   └── main.nf
│   ├── fgsea
│   │   └── main.nf
│   ├── salmon
│   │   └── main.nf
│   ├── trim_reads
│   │   └── main.nf
│   └── tximport
│       └── main.nf
├── nextflow.config
└── results

```
</details>

#### Overview
Here is a checklist and overview of important initial files in the project directory (and their relative paths) needed to run this demo along with descriptions of their origin:

`data/reference`:

- `metadata.csv`: contains information such as different IDs, sex, and treatments for the samples used for this demo. For this demo I selected 4 samples to subsample:

```
geo_accession,sex,treatment,SRR_ID,sample_type
GSM7265443,Female,Control,SRR24360653.sub,biol_rep_1
GSM7265449,Female,Control,SRR24360647.sub,biol_rep_2
GSM7265453,Female,Il10,SRR24360643.sub,biol_rep_1
GSM7265457,Female,Il10,SRR24360639.sub,biol_rep_2
```
- `cellchatv2_mouseLRI.rda`: CellChatDB v2 LRI annotations in an R Data File. When read in, it is a list of ligands and their corresponding list of receptors (in a list of lists format).

- `CellChatDB_preprocess.R`: **not covered in this demo** but in case there is interest, this is an Rscript details how the database was retrieved, cleaned up and transformed. 

`data/raw`:

- The following fastq files are subsampled versions of **raw paired-end `fastq` files** initially downloaded from the Sequence Read Archive (SRA) of the following [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA964442). These are archived  transcriptomics sequencing results for mouse samples exposed to one of two ligand moleculars to stimulate an immune response, a control group, and a treatment group that combines 1 ligand with a drug that blocks activity of the ligands, measured across 3 time points (n=20). For this purpose of this demo, I selected 2 samples from control and 2 exposed to only IL-10. See `metadata.csv` for further details.
- The original raw fastq were around 1.7 G per read file. These reads were subsampled as mentioned above using `seqtk`. For paired-end reads, it is recommended in the documentation to set a seed to ensure that the pairing is kept. I set seed = 100 and subsampled on 1 million reads.

<details> 
    <summary> Click here to see an example of how one fastq was generated using seqtk v1.4: </summary>

```bash
seqtk sample -s100 SRR24360639_1.fastq.gz 1000000 | gzip > SRR24360639.sub_1.fastq.gz
```
</details>
<br> 

`bin` is a directory that contains all Rscripts which the processes use:
- `tximport.R`
- `limma_voom.R`
- `plot_enhancedVolcano.R`   
- `ligand_enrichment_analysis.R`

`modules`: a directory of process definitions (placed into subdirectories) that are called into the main workflow script as detailed by the [Nextflow docs](https://www.nextflow.io/docs/latest/module.html).

`main.nf`: the main workflow script that will run the entire pipeline.

`results/`: the output directory where all processes instructed to write the outputs will go. The expected output are put into directories name after the process or tool that was used:
- `fastqc/`
- `ligand_enrichment/`
- `limma_voom/`
- `salmon_quant/`
- `trimmomatic/`
- `tximport/`
- `volcano_plot/`

## Quickstart
#### Run the pipeline:
To run the pipeline, make sure to be in the project directory. Then, run the `main.nf` as shown below.
```bash
cd ./biof501-term_project
nextflow run main.nf
```
This will run the entire workflow and a launch log will be displayed detailing the title of the workflow, some required inputs, and a list of processes it will run. Once a process is completed, it will show a checkmark to the right as shown below:

![Command Prompt Display](./assets/nextflow-startup.png)

### Step 1: Preprocess FASTQ files
### Step 2: Quantify reads using Salmon
### Step 3: Perform differential expression analysis
### Step 4: Perform pathway enrichment analysis
### Expected final results directory

## References
[1] Zhou L, Wang X, Peng L, et al. SEnSCA: Identifying possible ligand-receptor interactions and its application in cell-cell communication inference. J Cell Mol Med 28, e18372 (2024).

[2] Armingol E, Officer A, Harismendy O, et al. Deciphering cell–cell interactions and communication from gene expression. Nat Rev Genet 22, 71–88 (2021).

[3] Lim J, Park C, Kim M, et al. Advances in single-cell omics and multiomics for high-resolution molecular profiling. Exp Mol Med 56, 515–526 (2024).

[4] Sielemann K, Hafner A, Pucker B. The reuse of public datasets in the life sciences: potential risks and rewards. PeerJ. 8, e9954 (2020).

[5] Villemin JP, Bassaganyas L, Pourquier D, et al. Inferring ligand-receptor cellular networks from bulk and spatial transcriptomic datasets with BulkSignalR. Nucleic Acids Res 51, 4726–44 (2023). 

[6] Di Tommaso P, Chatzou M, Floden E, et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017).

[7] Ghosh, D. & Mersha, T. B. Publicly available cytokine data: Limitations and opportunities. J Allergy and Clin Immunol 150, 1053–1056 (2022). 

[8] Mootha VK, et al. Nat Genet 34, 267–273 (2003).

[9] Dimitrov D, Schäfer PSL, Farr E, et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. Nat Cell Biol 26, 1613–1622 (2024).

[10] Cui A, Huang T, Li S, et al. Dictionary of immune responses to cytokines at single-cell resolution. Nature 625, 377–384 (2024).

[11] Browaeys R, Saelens W, & Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods 17, 159–162 (2020).

[12] Jin, S., Plikus, M.V. & Nie, Q. CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nat Protoc (2024).

[13] fastqc

[14] trimmomatic

[15] salmon

[16] limma

[17] fgsea