include { FRED_PSEUDOUNIQ } from './modules/local/fred_remove_duplicates'
include { SAMTOOLS_FASTA } from './modules/local/samtools_fasta'
include { BLAST_INDEX_DB } from './modules/local/blast_indexdb'
include { BLAST_RUN } from './modules/local/blast_run'
include { BLAST2RMA } from './modules/local/blast2rma'
include { MEGAN_RMA2INFO } from './modules/local/megan_metagen'
include { SAMTOOLS_FQ2BAM } from './modules/local/samtools_import'
include { MAP_BWA } from './modules/local/map_bwa'
include { SAMTOOLS_FILTER } from './modules/local/samtools_filter'
include { SAMTOOLS_SORT } from './modules/local/samtools_sort'
include { BAM_RMDUP } from './modules/local/bam_rmdup'
include { BAM_DEAM_STATS } from './modules/local/bam_deam_stats'

// load the files

ch_split = Channel.fromPath("${params.split}/*.bam", checkIfExists:true)
ch_database = Channel.fromPath("${params.database}", checkIfExists:true)
ch_acc2taxid = Channel.fromPath("${params.acc2taxid}", checkIfExists:true)
ch_genomes = Channel.fromPath("${params.genomes}", checkIfExists:true)

ch_versions = Channel.empty()

workflow {

// add a fake meta
ch_split.map{it -> [['sample': it.baseName, 'id':it.baseName], it] }.set{ ch_split }


//
// 0. Remove Singletons (pseudouniq-step + filtering)
//


FRED_PSEUDOUNIQ(ch_split)

ch_pseudouniq = FRED_PSEUDOUNIQ.out.pseudouniq

ch_pseudouniq.map{meta, bam, stats -> 
    [
        meta+stats.splitCsv(sep:'\t', header:true).first(),
        bam
    ]
}.set{ ch_pseudouniq }

//
// 1. Convert bam to fasta (for blast)
//

SAMTOOLS_FASTA(ch_pseudouniq)

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

//
// 3. Run MEGAN
//

// here we need the fasta-file, the blastout and the acc2taxid-map
ch_blastout.combine(ch_fasta, by:0) //merge by meta
  .combine(ch_acc2taxid)
  .multiMap{meta, blastout, fasta, a2t -> 
    fasta: [meta, fasta]
    blastout: blastout
    a2t:a2t
  }
  .set{ch_b2rma_in}

// blast2rma

BLAST2RMA( ch_b2rma_in.fasta, ch_b2rma_in.blastout, ch_b2rma_in.a2t  )
ch_versions = ch_versions.mix(BLAST2RMA.out.versions.first())

ch_rma6 = BLAST2RMA.out.rma

// rma2info

MEGAN_RMA2INFO(ch_rma6)

ch_versions = ch_versions.mix(MEGAN_RMA2INFO.out.versions.first())

//
// NOW ITS SNAKEFILE p2
//

// We have the family-fasta-files in the 'fasta' channel

ch_family_fasta = MEGAN_RMA2INFO.out.fasta.transpose()
  .map{ meta, fasta ->
    [
      meta + ['family':fasta.baseName.split("-")[-1]],
      fasta
    ]
  }

// convert fasta to bam for mapping

SAMTOOLS_FQ2BAM(ch_family_fasta)
ch_versions = ch_versions.mix(SAMTOOLS_FQ2BAM.out.versions.first())

ch_bam = SAMTOOLS_FQ2BAM.out.bam.combine(ch_genomes)


// and map each family bam to all species in database under family/*.fasta

MAP_BWA(ch_bam)

// now filter MQ25
// run bam-rmdup
//--> goes into folder analyzed
// get substitution patterns
// print final report

MAP_BWA.out.mapped.transpose().map{ meta, bam ->
  def pattern = /^${meta.id}\.(.+)$/
  def matcher = bam.baseName =~ pattern
  [
    meta + ['species': matcher[0][1]],
    bam 
  ]
}.set{ch_mapped}

// Filter for Mapping quality
SAMTOOLS_FILTER ( ch_mapped )
ch_versions = ch_versions.mix( SAMTOOLS_FILTER.out.versions.first() )

ch_filtered = SAMTOOLS_FILTER.out.bam.map {
    [
      it[0]+['ReadsMapped':it[2].text.split(',')[1].trim() as int],
      it[1]
    ]
}

// Sort the bams
SAMTOOLS_SORT( ch_filtered )
ch_sorted = SAMTOOLS_SORT.out.bam

//
// BAM RM-DUP
//

BAM_RMDUP( ch_sorted )

// Add the Number of Unique Sequences to the meta
// first parse the bam-rmdup output file
stats = BAM_RMDUP.out.txt.splitCsv(sep:'\t', header:true, limit:1)

//now add the stats to the meta
BAM_RMDUP.out.bam
  .combine(stats, by: 0)
  .map{ meta, bam, stats ->
    [
        meta+["ReadsDeduped":stats["out"].replace(",","") as int],
        bam
    ]
  }
.set{ ch_deduped }

ch_versions = ch_versions.mix(BAM_RMDUP.out.versions.first())

BAM_DEAM_STATS( ch_deduped )

BAM_DEAM_STATS.out.tsv
  .map{ meta, tsv ->
      def deam = tsv.splitCsv(sep:'\t', header:true).first() // first because the splitCsv results have only one row -> [[key:value]]
      [
          meta+deam,
      ]
  }
.set{ ch_deam }


ch_deam.map{it[0]}.collectFile(
  name: "final_report.tsv", 
  newLine:true, 
  sort:true,
  storeDir:'.',
  seed:[
      "RG",
      "RawReads",
      "ReadsL35",
      "FinalNoPCRDups",
      "AvrgTimesSeen",
      "Family",
      "Species",
      "ReadsMapped",
      "ReadsDeduped",
      "Ancientness",
      "ReadsDeam1Term",
      "ReadsDeam3Term",
      "Deam5_95CI",
      "Deam3_95CI",
      "Deam5Cond_95CI",
      "Deam3Cond_95CI",
      "MeanFragmentLength",
      "MeanFragmentLength_Deam",
    ].join("\t"), 
  ){
    [
      it.id,
      it.all_seqs,
      it["seqs>=35"],
      it.final_noPCRdups_seqs,
      it.avrg_times_seen_final,
      it.family,
      it.species,
      it.ReadsMapped,
      it.ReadsDeduped,
      it.Ancientness,
      it['ReadsDeam(1term)'],
      it['ReadsDeam(3term)'],
      it['Deam5(95ci)'],
      it['Deam3(95ci)'],
      it['Deam5Cond(95ci)'],
      it['Deam3Cond(95ci)'],
      it['MeanFragmentLength'],
      it['MeanFragmentLength(3term)'],
    ].join('\t')
  }
}
