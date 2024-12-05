# BIOF 501 Term Project: RNASeq workflow for differential expression and pathway enrichment 

## Table of Content:
- [Background and Rationale](#background-and-rationale)
    - [Aims](#aims)
    - [Package Dependencies](#package-dependencies)
    - [Workflow Overview](#workflow-overview)
- [Usage](#usage)
    - [Installation](#installation)
    - [Step 1: Preprocess FASTQ files](#step-1-preprocess-fastq-files)
    - [Step 2: Quantify reads using Salmon](#step-2-quantify-reads)
    - [Step 3: Perform differential expression analysis](#step-3-perform-differential-expression-analysis)
    - [Step 4: Perform pathway enrichment analysis](#step-4-perform-pathway-enrichment-analysis)
- [References](#references)

## Background and Rationale
Cellular communication occurs through ligands binding to surface receptors, influencing signaling pathways that trigger phenotypic responses [[1](#references), [2](#references)]. Research on ligand-receptor interactions (LRIs) has gained momentum with single-cell omics, which provide high-resolution insights into individual cell states [[3](#references)]. However, bulk transcriptomic datasets, with their extensive clinical data from patient cohorts, remain valuable. Thus developing tools that complement single-cell studies using bulk transcriptomics to characterize cellular networks remain relevant [[4](#references)].

As omics technologies advance, complexity of data generated increase [[5](#references)]. In combination with complexity, rapid-evolving terminologies further create variability and barriers for reproducibility when leveraging these publicly available archival data sets[[6](#references)].  

### Aims:
**Aim 1 -** To provide a reproducible end-to-end workflow to analyze gene expression starting from raw fastq files.

**Aim 2 -** To explore gene set enrichment analysis in the context of LRIs characterized by pre-existing database annotations.

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
</details>

## Usage
### Base Installation Overview and Versions
The following versions of software were installed following the [Nextflow Docs]()
```
bash-5.2.21
docker-26.1.1
git-2.43.0
nextflow-24.10.0
open-jdk-17.0.10
```

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

[4] Villemin JP, Bassaganyas L, Pourquier D, et al. Inferring ligand-receptor cellular networks from bulk and spatial transcriptomic datasets with BulkSignalR. Nucleic Acids Res 51, 4726–44 (2023). 

[5] Di Tommaso P, Chatzou M, Floden E, et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017).

[6] Ghosh, D. & Mersha, T. B. Publicly available cytokine data: Limitations and opportunities. J Allergy and Clin Immunol 150, 1053–1056 (2022). 

## Initial Project directory structure
<details>
    <summary> Click here to see how the project directory should look like upon cloning. </summary>

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
