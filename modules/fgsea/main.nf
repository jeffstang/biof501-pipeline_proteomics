process PATHWAY_ENRICHMENT {
    tag "Perform Gene Set Enrichment on ${top_table_input.name} using ${database_of_interest.name}"
    publishDir "${params.outdir}/ligand_enrichment", mode: 'copy'

    input:
    path top_table_input
    path database_of_interest

    output:
    path "*.txt"

    script:
    """
    ligand_enrichment_analysis.R $top_table_input $database_of_interest
    """
}