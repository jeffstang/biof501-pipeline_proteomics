// This process downloads the fasta file
// based on the URL provided: demo fasta is from ENCODE GrCm39 transcriptome
// The reference files should be created in temp work directories

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
// provided some URL. The source should be the same as FASTA
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
