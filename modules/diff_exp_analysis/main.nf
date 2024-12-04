// This process runs the limma-voom package
process LIMMA_VOOM_DEA {
    publishDir "${params.outdir}/limma_voom", mode: "copy", overwrite: true
    input:
    path txi_input
    path metadata
    val  out_prefix

    output:
    path "*_topTable.csv", emit: logFC_topTable

    script:
    """
    limma_voom.R \
        $txi_input \
        $metadata \
        $out_prefix
    """
}