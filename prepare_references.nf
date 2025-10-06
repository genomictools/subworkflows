#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PFB }     from '../modules/bcftools/pfb.nf'
include { GCM }     from '../modules/penncnv/gcm.nf'
include { LEVELS }  from '../modules/quantisnp/levels.nf'

workflow prepare_references {
    take:
    dbsnp
    snplist
    gc
    tools

    main:
    dbsnp 
        | combine(snplist)
        | PFB
        | combine(gc)
        | GCM

    tools
        | filter { it == 'quantisnp' }
        | LEVELS

    emit:
    pfb = PFB.out
    gcm = GCM.out
    levels = LEVELS.out
}

workflow {
    dbsnp   = Channel.fromFilePairs(params.dbsnp, flat: true)
    snplist = Channel.fromPath(params.snplist)
    gc      = Channel.fromPath(params.gc)

    ref = prepare_references(dbsnp, snplist, gc)
}
