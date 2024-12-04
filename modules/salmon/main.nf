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