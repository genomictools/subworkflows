#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PFB } from '../modules/pfb.nf'
include { GCM } from '../modules/gcm.nf'

workflow prepare_references {
    take:
    dbsnp
    snplist
    gc
    
    main:
    dbsnp 
        | combine(snplist)
        | PFB
        | combine(gc)
        | GCM

    emit:
    pfb = PFB.out
    gcm = GCM.out
}

workflow {
    dbsnp   = Channel.fromFilePairs(params.dbsnp, flat: true)
    snplist = Channel.fromPath(params.snplist)
    gc      = Channel.fromPath(params.gc)

    ref = prepare_references(dbsnp, snplist, gc)
}
