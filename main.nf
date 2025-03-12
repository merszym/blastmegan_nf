include { BLAST_RUN } from './modules/local/blast_run'
include { BLAST_INDEX_DB } from './modules/local/blast_indexdb'
include { SAMTOOLS_FASTA } from './modules/local/samtools_fasta'


// load the files


ch_split = Channel.fromPath("${params.split}/*", checkIfExists:true)
ch_database = Channel.fromPath("${params.database}", checkIfExists:true)
ch_versions = Channel.empty()

workflow {

// add a fake meta
ch_split.map{it -> [['sample': it.baseName, 'id':it.baseName], it] }.set{ ch_split }

//
// 1. Convert bam to fasta
//

SAMTOOLS_FASTA(ch_split)

ch_fasta = SAMTOOLS_FASTA.out.fasta
ch_versions = ch_versions.mix(SAMTOOLS_FASTA.out.versions.first())

//
// 2. Run BLAST
//

// index the database
BLAST_INDEX_DB(ch_database)
ch_index = BLAST_INDEX_DB.out.index.toList()

// prepare channels for blast
ch_fasta.combine(ch_database)
    .multiMap{ meta, fasta, db ->
        fasta: [meta, fasta]
        database: db
    }
    .set{ch_blastin}

ch_blastdb = ch_blastin.database.combine(ch_index)

// Run Blast
BLAST_RUN(ch_blastin.fasta, ch_blastdb )

ch_blastout = BLAST_RUN.out.blastout
ch_versions = ch_versions.mix(BLAST_RUN.out.versions.first())

}