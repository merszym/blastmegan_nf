process FRED_MAKE_PSEUDOUNIQ{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label 'local'

    input:
    tuple val(meta), path(bam1), path(bam3), path(bam)

    output:
    tuple val(meta), path("masked_${bam1}"), path("masked_${bam3}"), path("masked_${bam}"), emit: bam

    script:
    """
    fred_remove_dups.py
    """
}