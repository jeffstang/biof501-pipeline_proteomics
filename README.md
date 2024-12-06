# BIOF 501 Term Project: RNASeq workflow for differential expression and pathway enrichment 

## Table of Contents:
- [Background and Rationale](#background-and-rationale)
    - [Objective and Aims](#objective-and-aims)
    - [Package Dependencies](#package-dependencies)
    - [Workflow Overview](#workflow-overview)
- [Usage](#usage)
    - [Installation](#installation)
    - [Quickstart](#quickstart)
    - [Project Structure and Overview](#project-structure-and-overview)
        - [Inputs](#overview-of-inputs)
        - [Important Expected Outputs](#overview-of-expected-important-outputs)
- [Detailed Description for Each Step](#detailed-description-for-each-step)
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
The workflow includes the following steps for analyzing mouse RNAseq samples:
1. Retrieve reference files from [GENCODE](https://www.gencodegenes.org/mouse/)
2. Preprocess raw fastq files
    - Provides quality metrics from `fastqc` [[13](#references)] pre- and post-trimming
    - Trims reads using `trimmomatic` [[14](#references)]
3. Count transcripts using `salmon` [[15](#references)] and generate a count matrix at the gene-level using `tximport` [[16](#references)]
4. Preprocess and analyze differentially expressed genes (DEGs) using `limma-voom` [[17](#references)] 
5. We leverage `fgsea` [[18](#references)] and an example database (CellChatDB v2 [[12](#references)]) to provide insight on ligand activity which can reflect on signalling pathway activity

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

### Quickstart
**Run the pipeline:**

To run the pipeline, make sure to be in the project directory. Then, run the `main.nf` as shown below.
```bash
cd ./biof501-term_project
nextflow run main.nf
```
This will run the entire workflow and a launch log will be displayed detailing the title of the workflow, some required inputs, and a list of processes it will run. Once a process is completed, it will show a checkmark to the right as shown below:

![Command Prompt Display](./assets/nextflow-startup.png)

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

#### Overview of Inputs
Here is a checklist and overview of important initial files in the project directory (and their relative paths) needed to run this demo along with descriptions of their origin:

`data/reference`:

- `metadata.csv`: contains information with fields including different IDs, sex, and treatments for the samples used for this demo (Control vs Il10 which is one of the ligands tested). For this demo I selected 4 samples to generate results for differential expression:

```
geo_accession,sex,treatment,SRR_ID,sample_type
GSM7265443,Female,Control,SRR24360653.sub,biol_rep_1
GSM7265449,Female,Control,SRR24360647.sub,biol_rep_2
GSM7265453,Female,Il10,SRR24360643.sub,biol_rep_1
GSM7265457,Female,Il10,SRR24360639.sub,biol_rep_2
```
- `cellchatv2_mouseLRI.rda`: CellChatDB v2 LRI annotations in an R Data File. When read in, it is a list of ligands and their corresponding list of receptors (in a list of lists format).

- `CellChatDB_preprocess.R`: **not covered in this demo** but in case there is interest, this is an Rscript detailing how the database was retrieved, cleaned up and transformed. 

`data/raw`:

- The following fastq files are subsampled versions of **raw paired-end `fastq` files** initially downloaded from the Sequence Read Archive (SRA) of the following [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA964442). These are archived  transcriptomics sequencing results for mouse samples exposed to one of two ligand moleculars to stimulate an immune response, a control group, and a treatment group that combines 1 ligand with a drug that blocks activity of the ligands, measured across 3 time points (n=20). For this purpose of this demo, I selected 2 samples from control and 2 exposed to only IL-10. See `metadata.csv` for further details.
- The original raw fastq were around 1.7 G per read file. These reads were subsampled as mentioned above using `seqtk`. For paired-end reads, it is recommended in the documentation to set a seed to ensure that the pairing is kept. I set seed = 100 and subsampled on 1 million reads.

<details> 
    <summary> Click here to see an example of how one fastq was generated using seqtk v1.4: </summary>

```bash
seqtk sample -s100 SRR24360639_1.fastq.gz 1000000 | gzip > SRR24360639.sub_1.fastq.gz
```
</details>

`bin/`: is a directory that contains all Rscripts named corresponding to the process used:
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

#### Overview of Expected Important Outputs
The following is an overview on the most important outputs to consider when looking at the results in the workflow. For a more detailed visual of outputs, please check the [Expected Final Results Directory](#expected-final-results-directory)

**FASTQC metrics:**

There should be two directories of FASTQC reports generated on each sample in `results/`. The first is named `raw/` and the second is named `trimmed/`. Each FASTQC sample directory should contain a .ZIP and .HTML file for R1 and R2. Opening an example .HTML file will show the full summary report for that fastq file like so:
![fastqc_report](./assets/fastqc_report.png)

**Differential Expression Analysis:**

A **top table** output as a "comma-separated value" (CSV) file that shows all genes and their significance scores ordered by most significant based on fitting the linear model from limma for treatment vs. control. Below is a subset example of the 5 most significant DEGs.

| Gene      | logFC   | AveExpr  | t          | P.Value           | adj.P.Val       | B          |
|-----------|---------|----------|------------|-------------------|-----------------|------------|
| Il4ra     | 5.3887  | 8.4107   | 26.6612    | 1.84e-09          | 1.47e-05        | 11.2194    |
| Ptpn1     | 1.8707  | 8.9710   | 19.1480    | 2.91e-08          | 1.16e-04        | 9.8171     |
| Plek      | 1.7910  | 8.7467   | 17.7162    | 5.53e-08          | 1.21e-04        | 9.1923     |
| Rab3il1   | 1.8201  | 8.4658   | 17.1868    | 7.10e-08          | 1.21e-04        | 8.9327     |
| Litaf     | 1.4867  | 9.2331   | 16.3916    | 1.05e-07          | 1.21e-04        | 8.5839     |

A **volcano plot** that provides a useful visual on genes that differ significantly between treatment compared to control.

![volcano_plot](./assets/volcano_plot.png)

**Ligand (Pathway) Enrichment Table:**

A table of significant ligands characterized by their receptor gene sets using CellChatDB LRI annotations ordered by their adjusted p values. Below shows a subset of the top 5 ligands represented. Note: Values may not make sense and have ties due to the subsetting of the fastq files. These may not make biological sense but it is expected that using the full raw fastq reads as input should generate meaningful results.

| Pathway  | pval           | padj           | log2err   | ES          | NES          | Size | Leading Edge |
|----------|----------------|----------------|-----------|-------------|--------------|------|--------------|
| Cd55     | 0.0015374361   | 0.0272538295   | 0.4551    | -0.9986     | -1.3357      | 1    | Adgre5       |
| Cd55b    | 0.0015374361   | 0.0272538295   | 0.4551    | -0.9986     | -1.3357      | 1    | Adgre5       |
| Cd86     | 0.0004179491   | 0.0272538295   | 0.4985    | -0.9997     | -1.3372      | 1    | Cd28         |
| Cntn1    | 0.0016150417   | 0.0272538295   | 0.4551    | -0.9877     | -1.4828      | 2    | Notch2       |
| Dlk1     | 0.0016150417   | 0.0272538295   | 0.4551    | -0.9877     | -1.4828      | 2    | Notch2       |

# Detailed Description For Each Step
## Step 1: Preprocess FASTQ files
Download Reference FASTA and GTF
  - Downloads the transcriptome reference assembly (sequence file) of GRCm39 from GENCODE and its basic annotation file that detailsdifferent IDs such as `gene_id` or`transcript_id`.  
  - The FASTA is used for indexing and quantifying transcripts 
  - The GTF is used for quantifying transcripts and converting transcript IDs into gene symbols

FastQC (Raw and Trimmed Reads):
- Analyze the quality of raw and trimmed sequencing reads
- Generate reports to detect issues like adapter contamination, GC bias,or low-quality reads [[13](#references)].
- This will ensure input reads meetquality standards for reliabletranscript counting

Trimmomatic
- Trims adapter sequences (which may have been introduced during library prep that are part of wet lab protocols with no biologicalrelevance) and low-quality bases from raw sequencing reads.
- For the `trimmomatic` tool, I set default flags and thresholds based on the [trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/TrimmomaticTrimmomaticManual_V0.32.pdf). Please adjust accordingly if needed.
- Improve the overall quality of reads for quantification.
- Expected output: are **paired trimmed reads** and **unpaired reads** for each mate 
    - Note that unpaired reads may be output if one mate (ie. R1) may have passed quality checks while the other read (R2) did not. Unpaired reads may still hold valuable information and can be used for single-read analyses or as additional data [[14](#references)] thus, I have output them into memory as well. 

The trimmed fastqs are the **input to be used for `salmon`**

### Step 2: Quantify reads using Salmon
Salmon (Index the FASTA)
- Salmon indexes the reference transcriptome (FASTA) in a format that is easy and fast to search. This step only needs to be done once and is not published into the results.
- Alternatively, you can also download pre-generated indices from other sources.

Salmon (Quantify Transcripts) [[15](#references)]
- Provides accurate and rapid transcript-level quantification essential for gene-level analysis
- The quantification step utilizes the indexed FASTA file and GTF annotations to map reads and sums up the paired trimmed fastqs per sample
- **IMPORTANT:** for this workflow, the relevant output that is used for downstream analyses are the **gene-level quantifications** named `quant.gene.sf` for each sample directory output. This is a text file with the following fields:
    - Name: transcript IDs 
    - Length: the length of the transcript in base pairs 
    - EffectiveLength: the computed effective length of the transcript
    - TPM: `salmon`'s estimates of the relative abundance of the transcript in transcripts per million
    - NumReads: number of reads mapped to each transcript.   

`tximport` [[16](#references)] is a tool that summarizes the above gene-level transcript counts per sample into a counts matrix format which is the preferred data format for downstream analyses like differential expression analysis. 
- The process `CREATE_TX2GENE` runs an `awk` command on the `GTF` annotation file that selects on the transcript ID and gene symbol to be output as a tab-delimited text file.
- This tx2gene file is used to sum up the transcripts at the gene-level and converts the names to human-readable gene symbols
- The final output is a gene x sample matrix where the rows represent all the genes captured for the demo data and the columns represent the SRR sample IDs. 

The output **gene-level counts matrix** of the `TXIMPORT` process is used for differential expression using `limma-voom` 

### Step 3: Perform differential expression analysis
Perform differential expression analysis using `limma`. Incorporates `voom` transformation to handle the mean-variance relationship of RNA-seq data [[17](#references)].

- The process takes in both the `tximport` generated counts matrix and the `metadata.csv` file
- The `metadata.csv` file contains the treatment required to set up the design matrix required for differential expression.

Expected output is a top table where the row names are the gene symbols, logFC represent the directionality of the gene expression (downregulated if negative or upregulated if positive), p-values represent significance, and the p-adjusted values are the significance values after multiple-test correction.

To easily visualize this explained output, a volcano plot was produced from this top table as shown [in the overview of important results.](#overview-of-expected-important-outputs)

### Step 4: Perform pathway enrichment analysis

Using the reference database of LRIs from CellChatDB v2, and the top table of DEGs from Step 3 as input and applying it to the `fgsea` tool [[18](#references)], I generated a ranked table of ligands enriched with differentially expressed receptor gene sets. An example of this table is also shown in the [overview of important results](#overview-of-expected-important-outputs). 

### Expected Final Results Directory

<details>
    <summary> Click here to see the tree structure of the expected results directory </summary>

```bash
results
├── fastqc
│   ├── raw
│   │   ├── SRR24360639.sub
│   │   │   ├── SRR24360639.sub_1_fastqc.html
│   │   │   ├── SRR24360639.sub_1_fastqc.zip
│   │   │   ├── SRR24360639.sub_2_fastqc.html
│   │   │   └── SRR24360639.sub_2_fastqc.zip
│   │   ├── ...
│   │   └── SRR24360653.sub
│   │       ├── SRR24360653.sub_1_fastqc.html
│   │       ├── SRR24360653.sub_1_fastqc.zip
│   │       ├── SRR24360653.sub_2_fastqc.html
│   │       └── SRR24360653.sub_2_fastqc.zip
│   └── trimmed
│       ├── SRR24360639.sub
│       │   ├── SRR24360639.sub_R1.trimmed_fastqc.html
│       │   ├── SRR24360639.sub_R1.trimmed_fastqc.zip
│       │   ├── SRR24360639.sub_R2.trimmed_fastqc.html
│       │   └── SRR24360639.sub_R2.trimmed_fastqc.zip
│       ├── ...
│       └── SRR24360653.sub
│           ├── SRR24360653.sub_R1.trimmed_fastqc.html
│           ├── SRR24360653.sub_R1.trimmed_fastqc.zip
│           ├── SRR24360653.sub_R2.trimmed_fastqc.html
│           └── SRR24360653.sub_R2.trimmed_fastqc.zip
├── ligand_enrichment
│   └── fgsea_ligand_activity.txt
├── limma_voom
│   └── logFC_DEG_topTable.csv
├── salmon_quant
│   ├── SRR24360639.sub
│   │   ├── aux_info
│   │   │   ├── ambig_info.tsv
│   │   │   ├── expected_bias.gz
│   │   │   ├── fld.gz
│   │   │   ├── meta_info.json
│   │   │   ├── observed_bias.gz
│   │   │   └── observed_bias_3p.gz
│   │   ├── cmd_info.json
│   │   ├── libParams
│   │   │   └── flenDist.txt
│   │   ├── lib_format_counts.json
│   │   ├── logs
│   │   │   └── salmon_quant.log
│   │   ├── quant.genes.sf
│   │   └── quant.sf
│   ├── ...
│   └── SRR24360653.sub
│       ├── aux_info
│       │   ├── ambig_info.tsv
│       │   ├── expected_bias.gz
│       │   ├── fld.gz
│       │   ├── meta_info.json
│       │   ├── observed_bias.gz
│       │   └── observed_bias_3p.gz
│       ├── cmd_info.json
│       ├── libParams
│       │   └── flenDist.txt
│       ├── lib_format_counts.json
│       ├── logs
│       │   └── salmon_quant.log
│       ├── quant.genes.sf
│       └── quant.sf
├── trimmomatic
│   ├── SRR24360639.sub
│   │   ├── SRR24360639.sub_R1.trimmed.fastq.gz
│   │   ├── SRR24360639.sub_R1.unpaired.fastq.gz
│   │   ├── SRR24360639.sub_R2.trimmed.fastq.gz
│   │   └── SRR24360639.sub_R2.unpaired.fastq.gz
│   └── SRR24360653.sub
│       ├── SRR24360653.sub_R1.trimmed.fastq.gz
│       ├── SRR24360653.sub_R1.unpaired.fastq.gz
│       ├── SRR24360653.sub_R2.trimmed.fastq.gz
│       └── SRR24360653.sub_R2.unpaired.fastq.gz
├── tximport
│   └── salmon_counts_txi.csv
└── volcano_plot
    └── volcano_plot.png
```
</details>

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

[16] Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4.

[17] limma

[18] fgsea