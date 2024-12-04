process ENHANCED_VOLCANO_PLOT {
    tag "Plotting Publication Ready Volcano Plot..."
    publishDir "${params.outdir}/volcano_plot", mode: 'copy'

    input:
    path toptable_input

    output:
    path "*.png"

    script:
    """
    plot_enhancedVolcano.R $toptable_input
    """
}