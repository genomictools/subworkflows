#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COMBINE }     from '../modules/gatk/combine.nf'
include { GENOTYPE }    from '../modules/gatk/genotype.nf'
include { MERGE }       from '../modules/gatk/merge.nf'
include { SPLIT }       from '../modules/bedtools/split.nf'

workflow joint_genotypes {
    take:
    cohorts

    main:
    // Load fasta
    fasta = Channel.fromPath(params.fasta) 
        | map { file -> [ params.assembly, file ] }
        | groupTuple(by: 0)

    // Create chunks from chrom sizes
    chunks =  Channel.of([params.assembly, file(params.chrom_sizes)])
        | SPLIT
        | transpose
        | map { assembly, bed -> [ assembly, "chunk_" + bed.name.split('\\.')[2], bed ] }

    // Collect by cohort and sample ID
    cohorts
        | combine(chunks, by: 0)
        | COMBINE
        | combine(fasta, by: 0)
        | GENOTYPE
        | groupTuple(by: [0, 1], sort: true)
        | combine(fasta, by: 0)
        | MERGE
        | set { joint_vcf_ch }

    emit:
    joint_vcf_ch
}

workflow {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.sample, row.file, row.index ]}

    joint_genotypes(cohorts_ch)
}
