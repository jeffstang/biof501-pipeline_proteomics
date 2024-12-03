/*
 * RNASeq Differential Analysis Pipeline
 * This is a pipeline designed for analyzing transcriptomics data in the context of cytokines. 
 * A fast approach of analyzing large RNAseq datasets starting from raw paired-end reads 
 * I prioritize computational efficiency using an alignment-free method to generate a report on cytokine perturbations.
 * 
 * Author: Jeffrey Tang
 *
 */

/*
 * Loading default parameters
*/

params.fastq = "$baseDir/data/raw/*{1,2}.fastq.gz"
params.fasta = "$baseDir/data/reference/grcm39_transcript.fa.gz"
params.gtf = "$baseDir/data/reference/grcm39_transcript.gtf.gz"
params.metadata_csv = "$baseDir/data/reference/metadata.csv" 
params.outdir = "results"

log.info """\
        RNASeq Differential Analysis Pipeline
        =====================================
        FASTQ                    : ${params.fastq}
        Transcriptome FASTA      : ${params.fasta}
        Gene Annotations         : ${params.gtf}
        Results Directory        : ${params.outdir}
        """
        .stripIndent()

fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.pc_transcripts.fa.gz"
gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.basic.annotation.gtf.gz"

//
// WORKFLOW: run main analysis pipeline
//

workflow {
    // Check if reference files (FASTA and GTF exist), download if necessary:
    DOWNLOAD_REFERENCES(fasta_url, gtf_url)
    
    // Initialize the read pair channel
    Channel
        .fromFilePairs( params.fastq, checkIfExists: true )
        .ifEmpty { error "Cannot find matching FASTQ files: ${params.fastq}" }
        .set { read_pairs_ch } 
         
    // Generate FASTQC reports on raw reads
    FASTQC( read_pairs_ch, "raw" )
 
    // Use trimmomatic to trim reads to remove low quality reads and adapter sequences
    TRIM_READS( read_pairs_ch )
        
    // Run Salmon processes
    SALMON_INDEX( DOWNLOAD_REFERENCES.out.fasta )
    SALMON_QUANT( 
        SALMON_INDEX.out.idx, 
        DOWNLOAD_REFERENCES.out.annotation, 
        TRIM_READS.out.trimmed_paired
        )

    // Create a tx2gene file from annotations GTF 
    CREATE_TX2GENE(DOWNLOAD_REFERENCES.out.annotation)

    // Gather all salmon-processed samples and define necessary inputs
    quant_dirs = SALMON_QUANT.out.quant.collect()
    tx2gene_input = CREATE_TX2GENE.out.tx2gene
    
    // Use tximport to merge and convert transcript IDs to gene names
    // This will generate a gene by samples counts matrix
    TXIMPORT_PROCESS(
    quant_dirs,
    tx2gene_input,
    'salmon_counts'
    )

    // Run differential expression analysis using limma-voom from tximport_process
    // output, which is a single salmon counts matrix. 
    // The inputs for this process are the tximport object, the metadata, and the output prefix to name the CSV
    // The expected output is the top table of differentially expressed genes
    txi_input = TXIMPORT_PROCESS.out.txi_object
    
    LIMMA_VOOM_DEA( 
        txi_input,
        file(params.metadata_csv),
        "logFC_DEG"
    )

}

/* 
 * Subworkflows:
 * Download FASTA and GTF files
*/
workflow DOWNLOAD_REFERENCES {
    take:
    fasta_url
    gtf_url
    
    main:
    // Retrieve transcript fasta
    DOWNLOAD_FASTA(fasta_url)

    // Retrieve transcript annotations file
    DOWNLOAD_GTF(gtf_url)

    emit:
    fasta = DOWNLOAD_FASTA.out.fasta
    annotation = DOWNLOAD_GTF.out.annotation
}

// Processes:
// This process downloads the fasta file
process DOWNLOAD_FASTA {
    errorStrategy 'retry' 
    maxRetries 2
    cache 'deep'
    
    input:
    val(fasta_url)

    output: 
    path "*.fa.gz", emit: fasta
                
    script:
    """
    wget -O ${fasta_url.split("/")[-1]} ${fasta_url}
    """    
}

// This process downloads the GTF annotation file
process DOWNLOAD_GTF {
    errorStrategy 'retry'
    maxRetries 2
    cache 'deep'

    input:
    val(gtf_url)

    output: 
    path "*.gtf.gz", emit: annotation

    script:
    """
    wget -O ${gtf_url.split("/")[-1]} ${gtf_url}
    """
}

// This process reads the fastq files and outputs the QC metrics from the fastqc tool
process FASTQC {    
    tag "FASTQC on ${read_type} reads for sample ${sample_id}"
    publishDir "${params.outdir}/fastqc/${read_type}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)    
    val(read_type)
    
    output:
    path "*_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    fastqc -f fastq -q ${reads} 
    """
}

// This is a cloned process of FASTQC that is used specifically to report 
// the QC metrics after trimming the reads
process TRIMMED_FASTQC {    
    tag "FASTQC on ${read_type} reads for sample ${sample_id}"
    publishDir "${params.outdir}/fastqc/${read_type}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)    
    val(read_type)
    
    output:
    path "*_fastqc.{zip,html}", emit: trimmed_fastqc_reports

    script:
    """
    fastqc -f fastq -q ${reads} 
    """
}

// Uses trimmomatic to trim the reads, it will process trimmed fastqs 
// and a separate file for unpaired reads
process TRIM_READS {
    tag "TRIM_READS on ${sample_id}"
    publishDir "${params.outdir}/trimmomatic/${sample_id}", mode: 'copy'
    
    input: 
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz"), emit: trimmed_paired
    tuple val(sample_id), path("${sample_id}_R{1,2}.unpaired.fastq.gz"), emit: trimmed_unpaired
    
    script:
    def (r1, r2) = reads
    """
    trimmomatic PE -threads 4 \
        $r1 $r2 \
        ${sample_id}_R1.trimmed.fastq.gz \
        ${sample_id}_R1.unpaired.fastq.gz \
        ${sample_id}_R2.trimmed.fastq.gz \
        ${sample_id}_R2.unpaired.fastq.gz \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Indexes the fasta file so that it can be used for the salmon tool
process SALMON_INDEX {
    tag "Creating Salmon index from ${transcriptome_fasta.name}"
    input:
    path(transcriptome_fasta)
    
    output:
    path "salmon_index", emit: idx

    script:
    """
    salmon index --threads 4 -t $transcriptome_fasta -i salmon_index --gencode
    """
}

// Processes the trimmed fastq files and sums up transcript and gene-level counts
process SALMON_QUANT {
    tag "Converting counts from ${sample_id}"
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    path index
    path gtf_gz
    tuple val(sample_id), path(reads)

    output:
    path(sample_id), emit: quant

    script:
    def (r1, r2) = reads
    """
    gunzip -c ${gtf_gz} > annotation.gtf
    salmon quant \
    --threads 4 --libType A \
    --index $index \
    --validateMappings \
    --geneMap annotation.gtf \
    --output ${sample_id} \
    -1 ${r1} \
    -2 ${r2}
    rm annotation.gtf
    """
}

// Create a transcript to gene annotation file that's tab-delimited
// Used to convert transcript IDs to gene symbols
process CREATE_TX2GENE {
    input:
    path gtf_gz

    output:
    path "tx2gene.tsv", emit: tx2gene
    
    script:
    """
    gunzip -c ${gtf_gz} | awk '\$3 == "gene" { 
        match(\$0, /gene_id "[^"]+"/); 
        gene_id = substr(\$0, RSTART+9, RLENGTH-10); 
        match(\$0, /gene_name "[^"]+"/); 
        gene_name = substr(\$0, RSTART+11, RLENGTH-12); 
        print gene_id, gene_name; 
    }' > tx2gene.tsv
    """
}

// Process salmon quant files and merges samples into a single matrix
// This creates a txi object that contains abundance, counts, and transcript length information
// Although limma-voom does not require $length, I added this for modularity in case
// I wish to run DeSeq2 in the future
process TXIMPORT_PROCESS {
    publishDir "${params.outdir}/tximport", mode: 'copy', overwrite: true

    input:
    path quantified_gene_counts
    path tx2gene_tsv
    val  output_prefix

    output:
    path "${output_prefix}_txi.csv", emit: txi_object
    script:
    """
    tximport.R \
        "$quantified_gene_counts" \
        $tx2gene_tsv \
        $output_prefix
    """
}

// This process runs the limma-voom package
process LIMMA_VOOM_DEA {
    publishDir "${params.outdir}/limma_voom", mode: "copy", overwrite: true
    input:
    path txi_input
    path metadata
    val  out_prefix

    output:
    path "${output_prefix}_topTable.csv", emit: logFC_topTable

    script:
    """
    limma_voom.R \
        $txi_input \
        $metadata \
        $out_prefix
    """
}