process MAP_BWA {
    container (workflow.containerEngine ? "merszym/network-aware-bwa:v0.5.10" : null)
    tag "${meta.id}:${meta.family}"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(bam), path(genome)

    output:
    tuple val(meta), path("*.bam"), emit: mapped
    path 'versions.yml'           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """

    for fasta in \$(ls ${genome}/${meta.family}/*.fasta); do
        filename="\${fasta##*/}"
        speciesname="\${filename%.fasta}"
        bwa bam2bam -g "\$fasta" ${args} --only-aligned ${bam} > "${meta.sample}.\$speciesname".bam;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 > /dev/null | grep Version | cut -d ' ' -f2)
    END_VERSIONS
    """
}