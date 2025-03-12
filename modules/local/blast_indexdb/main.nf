process BLAST_INDEX_DB{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.9.0--pl526he19e7b1_7' :
        'quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7' }"

    input:
    path('database.fasta')

    output:
    path("database.fasta.*"), emit: index

    script:
    """
    makeblastdb -in database.fasta -dbtype nucl -parse_seqids
    """
}