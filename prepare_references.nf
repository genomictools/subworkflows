#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PFB }         from '../modules/bcftools/pfb.nf'
include { GCM }         from '../modules/penncnv/gcm.nf'

workflow prepare_references {
    take:
    snplist
    dbsnp
    gc

    main:
    snplist = Channel.fromPath(params.snplist)
    dbsnp   = Channel.fromFilePairs(params.dbsnp, flat: true)
    gc      = Channel.fromPath(params.gc)

    snplist 
        | splitText(by: params.chunk, file: true)
        | map { file -> tuple( "snps_${file.name.tokenize('\\.')[-2].toInteger()}", file ) }
        | combine(dbsnp)
        | PFB
        | combine(gc)
        | GCM

    PFB.out | collectFile(keepHeader: true) | set { pfb }
    GCM.out | collectFile(keepHeader: true) | set { gcm }

    emit:
    pfb
    gcm
}

workflow {
    ref = prepare_references(params.snplist, params.dbsnp, params.gc)
}
