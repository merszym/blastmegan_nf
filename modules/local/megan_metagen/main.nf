process MEGAN_RMA2INFO{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megan:6.25.10--h9ee0642_0' :
        'quay.io/biocontainers/megan:6.25.10--h9ee0642_0' }"
    tag "$meta.id"

    input:
    tuple val(meta), path(rma6)

    output:
    tuple val(meta), path("${meta.sample}.tsv"), emit: tsv
    tuple val(meta), path("*.fasta")           , emit: fasta 
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    rma2info -i ${rma6} -l -c2c Taxonomy -r -n -mro -s \
      |  grep -E "^F" | cut -f 2,3 > ${meta.sample}.tsv
    

    for fam in \$(cut -f1 ${meta.sample}.tsv); do
        read-extractor -i ${rma6} -c Taxonomy -n \$fam --allBelow > ${meta.sample}_reads-\$fam.fasta;
    done
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rma2info/read-extractor: \$(echo "MEGAN CE (v 6.25.10)")
    END_VERSIONS
    """
}