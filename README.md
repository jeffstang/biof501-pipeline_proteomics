# BIOF 501 Term Project: RNASeq workflow for differential expression and gene summarization 

## Background and Rationale
• include the what’s and why’s – also your aims

• include any package dependencies that are required (bullet points are ok for this)

• You can include your DAG here

## Usage
Make sure you format everything so that step by step usage details are included. ## If we can’t run your pipeline then we can’t give you marks.

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