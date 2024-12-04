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