#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COMBINE }     from '../modules/gatk/combine.nf'
include { GENOTYPE }    from '../modules/gatk/genotype.nf'
include { MERGE }       from '../modules/gatk/merge.nf'
include { SPLIT }       from '../modules/bedtools/split.nf'

workflow joint_genotypes {
    take:
    samplesheet
    gvcf

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
        | take(3)

    // Collect by cohort and sample ID
    gvcf
        | combine(samplesheet, by: 0)
        | groupTuple(by: [0,2])
        | filter { id, files, cohort -> files.flatten().size() == 2 }
        | map { id, files, cohort -> [ params.assembly, cohort, id ] + files.flatten().sort { it.name } }
        | take(3)
        | groupTuple(by: [0, 1])
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
    samplesheet_ch = Channel.fromPath(params.samplesheet)
        | splitCsv(header: false, sep: '\t', skip: 1)
        | map { row -> [ params.cohort, row.name] }

    downloaded_files = Channel.fromPath(params.files)
        | splitCsv(header: false, sep: '\t', skip: 1)
        | map { row -> [ row.id, row.file] }

    joint_genotypes(samplesheet_ch, downloaded_files)
}
