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