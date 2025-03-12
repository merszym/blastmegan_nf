process BLAST_RUN{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.9.0--pl526he19e7b1_7' :
        'quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7' }"
    tag "$meta.id"
    publishDir "blast", mode: 'copy', pattern: '*.blastout'

    input:
    tuple val(meta), path(fasta)
    tuple path('database.fasta'), path(index)

    output:
    tuple val(meta), path("${meta.sample}.blastout"), emit: blastout
    path "versions.yml"                             , emit: versions

    script:
    """
    blastn -query ${fasta} -db database.fasta -outfmt 6 -num_threads ${params.threads} -out ${meta.sample}.blastout

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -h | tail -3 | head -1)
    END_VERSIONS
    """
}