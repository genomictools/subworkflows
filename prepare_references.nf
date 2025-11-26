#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CAF }         from '../modules/bcftools/caf.nf'
include { PFB }         from '../modules/rocker/pfb.nf'
include { GCM }         from '../modules/penncnv/gcm.nf'

workflow prepare_references {
    take:
    snplist
    dbsnp
    gc

    main:
    snplist = Channel.fromPath(params.snplist)
        | splitText(by: params.chunk, file: true)
        | map { file -> tuple( "snps_${file.name.tokenize('\\.')[-2].toInteger()}", file ) }
        
    dbsnp   = Channel.fromFilePairs(params.dbsnp, flat: true)
        | CAF
    gc      = Channel.fromPath(params.gc)

    snplist 
        | combine(dbsnp)
        | PFB
        | groupTuple(by: 0)
        | combine(gc)
        | GCM

    PFB.out | collectFile(keepHeader: true, storeDir: "${params.output_dir}/ref") { [ "${it[0]}.pfb", it[2] ] } | set { pfb }
    GCM.out | collectFile(keepHeader: true, storeDir: "${params.output_dir}/ref") { [ "${it[0]}.gcm", it[2] ] } | set { gcm }

    emit:
    pfb
    gcm
}

workflow {
    ref = prepare_references(params.snplist, params.dbsnp, params.gc)
}
