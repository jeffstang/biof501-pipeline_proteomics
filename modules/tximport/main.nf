// Process salmon quant files and aggregates sample counts into a single counts matrix
// Of note: I output a CSV for the purposes of simplicity for the project
// Future improvements can be made such that it creates and RDS that stores the TXI object

process TXIMPORT_PROCESS {
    publishDir "${params.outdir}/tximport", mode: 'copy', overwrite: true

    input:
    // inputs include:
    // 1) the salmon_quant result directory, 
    // 2) the transcript to gene file derived from the GTF 
    // 3) an output prefix to name the output file
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