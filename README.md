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
Cellular communication is dictated by the release of molecules called ligands that bind to cellular surface receptors to illicit a response [1]. These ligand-receptor interactions (LRIs) can be deterministic of signalling pathway activity [2].

### Aims
Aim 1: To provide a generic 

### Package Dependencies

### Workflow Overview
The workflow includes the following steps:

1. Download Reference files (FASTA and GTF)
2. Preprocess raw fastq files (includes generating QC reports via `fastqc` and trimming the reads using `trimmomatic`)
3. Count transcripts and generate count matrix using `salmon`
4. Analyze differentially expressed genes (DGE) using `limma-voom`
5. Provide insight on potential pathways relevant to significant DGEs using `fgsea`

<details>
    <summary>You can include find the DAG workflow here</summary>
![Workflow DAG](path/to/DAG-image.png)

## Usage
### Installation
### 1.
### 2.
### 3.

• Make sure you format everything so that step by step usage details are included. 

## If we can’t run your pipeline then we can’t give you marks.

• Installation (if necessary) including any datasets that are to be used if they are not provided (i.e. how to download them using wget or curl – exact paths need to be specified and the data must be accessible)

• Exact step by step usage with descriptive comments on what action is being performed in each step

## Project directory structure
<details>
    <summary> Click here to see the drop-down view of project directory in a tree-like format. </summary>
```bash
├── README.md
├── bin
│   ├── limma_voom.R
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
│       └── metadata.csv
├── main.nf
├── modules
│   ├── convert_tx2gene
│   │   └── main.nf
│   ├── diff_exp_analysis
│   │   └── main.nf
│   ├── download_ref
│   │   └── main.nf
│   ├── fastqc
│   │   └── main.nf
│   ├── salmon
│   │   └── main.nf
│   ├── trim_reads
│   │   └── main.nf
│   └── tximport
│       └── main.nf
├── nextflow.config
├── results
│   ├── fastqc
│   │   └── raw
│   │       ├── SRR24360639.sub
│   │       │   ├── SRR24360639.sub_1_fastqc.html
│   │       │   ├── SRR24360639.sub_1_fastqc.zip
│   │       │   ├── SRR24360639.sub_2_fastqc.html
│   │       │   └── SRR24360639.sub_2_fastqc.zip
│   │       ├── SRR24360643.sub
│   │       │   ├── SRR24360643.sub_1_fastqc.html
│   │       │   ├── SRR24360643.sub_1_fastqc.zip
│   │       │   ├── SRR24360643.sub_2_fastqc.html
│   │       │   └── SRR24360643.sub_2_fastqc.zip
│   │       ├── SRR24360647.sub
│   │       │   ├── SRR24360647.sub_1_fastqc.html
│   │       │   ├── SRR24360647.sub_1_fastqc.zip
│   │       │   ├── SRR24360647.sub_2_fastqc.html
│   │       │   └── SRR24360647.sub_2_fastqc.zip
│   │       └── SRR24360653.sub
│   │           ├── SRR24360653.sub_1_fastqc.html
│   │           ├── SRR24360653.sub_1_fastqc.zip
│   │           ├── SRR24360653.sub_2_fastqc.html
│   │           └── SRR24360653.sub_2_fastqc.zip
│   ├── limma_voom
│   │   └── logFC_DEG_topTable.csv
│   ├── salmon_quant
│   │   ├── SRR24360639.sub
│   │   │   ├── aux_info
│   │   │   │   ├── ambig_info.tsv
│   │   │   │   ├── expected_bias.gz
│   │   │   │   ├── fld.gz
│   │   │   │   ├── meta_info.json
│   │   │   │   ├── observed_bias.gz
│   │   │   │   └── observed_bias_3p.gz
│   │   │   ├── cmd_info.json
│   │   │   ├── libParams
│   │   │   │   └── flenDist.txt
│   │   │   ├── lib_format_counts.json
│   │   │   ├── logs
│   │   │   │   └── salmon_quant.log
│   │   │   ├── quant.genes.sf
│   │   │   └── quant.sf
│   │   ├── SRR24360643.sub
│   │   │   ├── aux_info
│   │   │   │   ├── ambig_info.tsv
│   │   │   │   ├── expected_bias.gz
│   │   │   │   ├── fld.gz
│   │   │   │   ├── meta_info.json
│   │   │   │   ├── observed_bias.gz
│   │   │   │   └── observed_bias_3p.gz
│   │   │   ├── cmd_info.json
│   │   │   ├── libParams
│   │   │   │   └── flenDist.txt
│   │   │   ├── lib_format_counts.json
│   │   │   ├── logs
│   │   │   │   └── salmon_quant.log
│   │   │   ├── quant.genes.sf
│   │   │   └── quant.sf
│   │   ├── SRR24360647.sub
│   │   │   ├── aux_info
│   │   │   │   ├── ambig_info.tsv
│   │   │   │   ├── expected_bias.gz
│   │   │   │   ├── fld.gz
│   │   │   │   ├── meta_info.json
│   │   │   │   ├── observed_bias.gz
│   │   │   │   └── observed_bias_3p.gz
│   │   │   ├── cmd_info.json
│   │   │   ├── libParams
│   │   │   │   └── flenDist.txt
│   │   │   ├── lib_format_counts.json
│   │   │   ├── logs
│   │   │   │   └── salmon_quant.log
│   │   │   ├── quant.genes.sf
│   │   │   └── quant.sf
│   │   └── SRR24360653.sub
│   │       ├── aux_info
│   │       │   ├── ambig_info.tsv
│   │       │   ├── expected_bias.gz
│   │       │   ├── fld.gz
│   │       │   ├── meta_info.json
│   │       │   ├── observed_bias.gz
│   │       │   └── observed_bias_3p.gz
│   │       ├── cmd_info.json
│   │       ├── libParams
│   │       │   └── flenDist.txt
│   │       ├── lib_format_counts.json
│   │       ├── logs
│   │       │   └── salmon_quant.log
│   │       ├── quant.genes.sf
│   │       └── quant.sf
│   ├── trimmomatic
│   │   ├── SRR24360639.sub
│   │   │   ├── SRR24360639.sub_R1.trimmed.fastq.gz
│   │   │   ├── SRR24360639.sub_R1.unpaired.fastq.gz
│   │   │   ├── SRR24360639.sub_R2.trimmed.fastq.gz
│   │   │   └── SRR24360639.sub_R2.unpaired.fastq.gz
│   │   ├── SRR24360643.sub
│   │   │   ├── SRR24360643.sub_R1.trimmed.fastq.gz
│   │   │   ├── SRR24360643.sub_R1.unpaired.fastq.gz
│   │   │   ├── SRR24360643.sub_R2.trimmed.fastq.gz
│   │   │   └── SRR24360643.sub_R2.unpaired.fastq.gz
│   │   ├── SRR24360647.sub
│   │   │   ├── SRR24360647.sub_R1.trimmed.fastq.gz
│   │   │   ├── SRR24360647.sub_R1.unpaired.fastq.gz
│   │   │   ├── SRR24360647.sub_R2.trimmed.fastq.gz
│   │   │   └── SRR24360647.sub_R2.unpaired.fastq.gz
│   │   └── SRR24360653.sub
│   │       ├── SRR24360653.sub_R1.trimmed.fastq.gz
│   │       ├── SRR24360653.sub_R1.unpaired.fastq.gz
│   │       ├── SRR24360653.sub_R2.trimmed.fastq.gz
│   │       └── SRR24360653.sub_R2.unpaired.fastq.gz
│   └── tximport
│       └── salmon_counts_txi.csv
└── run.sh
```

## Reference
1. 