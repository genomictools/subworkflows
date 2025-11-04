#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COMBINE }     from '../modules/gatk/combine.nf'
include { GENOTYPE }    from '../modules/gatk/genotype.nf'
include { MERGE }       from '../modules/gatk/merge.nf'
include { createChunks }from '../modules/utils/creatchunks.nf'

workflow joint_genotypes {
    take:
    samplesheet
    gvcf

    main:
    // Load fasta
    fasta = Channel.fromPath(params.fasta) 
        | map { file -> [ file.simpleName, file ] }
        | groupTuple(by: 0)

    // Create chunks from chrom sizes
    chunks =  createChunks(params.chrom_sizes, params.chunk) | take(3)

    // Collect by cohort and sample ID
    gvcf
        | groupTuple(by: 0, sort: true)
        | combine(samplesheet, by: 0)
        | map { id, files, cohort -> [ cohort, id ] + files.flatten() }
        | groupTuple(by: 0)
        | combine(fasta)
        | combine(chunks)
        | COMBINE
        | combine(fasta)
        | GENOTYPE
        | groupTuple(by: [0, 1], sort: true)
        | MERGE 
        | set { joint_vcf_ch }

    emit:
    joint_vcf_ch
}

workflow {
    joint_genotypes(samplesheet_ch, gvcf_ch, fasta_ch)
}
