process BLAST2RMA{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megan:6.25.10--h9ee0642_0' :
        'quay.io/biocontainers/megan:6.25.10--h9ee0642_0' }"
    
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)
    path(blastout)
    path(accession2taxid)

    output:
    tuple val(meta), path("${meta.sample}.rma6"), emit: rma
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    blast2rma \
       -f BlastText \
       -bm BlastN \
       -i ${blastout} \
       -r ${fasta} \
       -a2t ${accession2taxid} \
       $args \
       -o ${meta.sample}.rma6

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast2rma: \$(echo "MEGAN CE (v 6.25.10)")
    END_VERSIONS
    """
}