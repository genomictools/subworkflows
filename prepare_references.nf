#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PFB }         from '../modules/bcftools/pfb.nf'
include { GCM }         from '../modules/penncnv/gcm.nf'
include { QUANTMODEL }  from '../modules/rocker/quantmodel.nf'

workflow prepare_references {
    take:
    snplist
    dbsnp
    gc

    main:
    // Initialize empty channels
    pfb = gcm = quantlevels = quantparams = Channel.empty()

    if ( params.tools.contains('penncnv') ) {
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
    }
    
    if ( params.tools.contains('quantisnp') ) {
        Channel.of('levels', 'params')
            | QUANTMODEL
            | branch {
                quantlevels: it.first() == 'levels'
                quantparams: it.first() == 'params'
            }
            | set { models }
        models.quantlevels | collectFile(keepHeader: true) | set { quantlevels }
        models.quantparams | collectFile(keepHeader: true) | set { quantparams }
    }
    emit:
    pfb
    gcm
    quantlevels
    quantparams
}

workflow {
    ref = prepare_references(params.snplist, params.dbsnp, params.gc)
}
