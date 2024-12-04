// Create a transcript to gene annotation file that's tab-delimited
// Used to convert gene IDs to gene symbols
// The purpose of this is to translate gene symbols to human readable gene symbols
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