#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CAF } from '../modules/bcftools/caf.nf'
include { PFB } from '../modules/rocker/pfb.nf'
include { GCM } from '../modules/penncnv/gcm.nf'

workflow prepare_references {
    take:
    snplist
    known_sites
    dbsnp
    gc

    main:
    // Extract CAF info from dbsnp
    Channel.fromFilePairs(dbsnp, flat: true)
        | CAF
        | set { dbsnp }

    // Load GC content file
    Channel.fromPath(gc) | set { gc }

    // Prepare PFB file
    Channel.empty()
        | concat(Channel.fromPath(snplist))
        | ( params.known_sites != null ? concat(Channel.fromPath(known_sites)) : identity() )
        | splitText(by: params.chunk, file: true)
        | map { file -> tuple( "${file.name.tokenize('\\.')[0]}.${file.name.tokenize('\\.')[-2].toInteger()}", file ) }
        | combine(dbsnp)
        | PFB
        | collectFile(
            keepHeader: true,
            storeDir: "${params.output_dir}/ref"
        ) { [ "${it[0]}.pfb", it[2] ] }
        | set { pfb }

    // Prepare GCM file
    PFB.out
        | groupTuple(by: 0)
        | combine(gc)
        | GCM
        | collectFile(
            keepHeader: true,
            storeDir: "${params.output_dir}/ref"
        ) { [ "${it[0]}.gcm", it[2] ] }
        | set { gcm }

    emit:
    pfb
    gcm
}

workflow {
    ref = prepare_references(params.snplist, params.known_sites, params.dbsnp, params.gc)
}
