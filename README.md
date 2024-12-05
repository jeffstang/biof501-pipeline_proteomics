# BIOF 501 Term Project: RNASeq workflow for differential expression and pathway enrichment 

## Table of Contents:
- [Background and Rationale](#background-and-rationale)
    - [Objective and Aims](#objective-and-aims)
    - [Package Dependencies](#package-dependencies)
    - [Workflow Overview](#workflow-overview)
- [Usage](#usage)
    - [Installation](#installation-and-overview)
    - [Initial Project Directory Structure](#initial-project-directory-structure)
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

**Aim 2 -** To explore gene set enrichment analysis in the context of LRIs described by pre-existing database annotations. A **key assumption** of note is that we receptor gene exprression changes can characterize the ligand and sufficiently describes the "pathway" of ligand activity.  

### Workflow Overview
The workflow includes the following steps:

1. Retrieve reference files from ENCODE
2. Preprocess raw fastq files  (provides quality metrics from `fastqc` pre- and post-trimming and trims reads using `trimmomatic`)
3. Count transcripts and generate count matrix using `salmon`
4. Preprocess and analyze differentially expressed genes (DGE) using `limma-voom` 
5. We leverage `fgsea` and the example database to provide insight on ligand activity which can reflect on signalling pathway activity

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
        DEA --> VOLCANO[ENHANCED_VOLCANO_PLOT: Generate a volcano plot to show DGEs]
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
#### List of docker images used (see[ nextflow.config](https://github.com/jeffstang/biof501-term_project/blob/main/nextflow.config)):
```
edgeR
quay.io/biocontainers/bioconductor-edger:4.0.2--r43hf17093f_0

Enhanced Volcano
quay.io/biocontainers/bioconductor-enhancedvolcano:1.20.0--r43hdfd78af_0

FastQC
staphb/fastqc:0.12.1

fgsea
quay.io/biocontainers/bioconductor-fgsea:1.28.0--r43hf17093f_1 

salmon
combinelab/salmon:1.10.3

trimmomatic
quay.io/biocontainers/trimmomatic:0.36--4                 

tximport
quay.io/biocontainers/bioconductor-tximport:1.26.0--r42hdfd78af_0
```

## Usage
### Installation
The following versions of software were installed following the [Nextflow Docs]()
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

### Initial Project Directory Structure
Upon cloning, the repository should come with a set of base files.
<details>
    <summary> Click here to see a tree overview of the overall directory structure </summary>

```bash
tree biof501-term_project
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
├── results
└── run.sh
```
</details>

Here is a checklist and details on important files:

- `data/reference`
    - `metadata.csv`: contains information such as different IDs, sex, and treatments for the samples used for this tutorial
    - `cellchatv2_mouseLRI.rda`: 
- 

### Step 1: Preprocess FASTQ files
### Step 2: Quantify reads using Salmon
### Step 3: Perform differential expression analysis
### Step 4: Perform pathway enrichment analysis
### Expected final results directory

• Installation (if necessary) including any datasets that are to be used if they are not provided (i.e. how to download them using wget or curl – exact paths need to be specified and the data must be accessible)

• Exact step by step usage with descriptive comments on what action is being performed in each step

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