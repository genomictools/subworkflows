#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COMBINE }     from '../modules/gatk/combine.nf'
include { GENOTYPE }    from '../modules/gatk/genotype.nf'
include { GATHER }      from '../modules/gatk/gather.nf'

workflow joint_genotypes {
    take:
    samplesheet
    gvcf
    fasta

    main:
    gvcf
        | groupTuple(by: 0, sort: true)
        | combine(samplesheet, by: 0)
        | map { id, files, cohort -> [ cohort, id ] + files.flatten() }
        | groupTuple(by: 0)
        | combine(fasta)
        | COMBINE // Make optional
        | combine(fasta, by: [0, 1])
        | GENOTYPE // Make optional
        | groupTuple(by: [0,2], sort: true)
        | GATHER // Make optional
        // | set { joint_vcf_ch }

    // emit:
    // joint_vcf_ch
    | view
}

workflow {
    joint_genotypes(samplesheet_ch, gvcf_ch, fasta_ch)
}
