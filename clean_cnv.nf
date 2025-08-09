#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FILTER }    from '../modules/filter.nf'
include { CLEAN }     from '../modules/clean.nf'
include { SCAN }      from '../modules/scan.nf'
include { VISUALIZE } from '../modules/visualize.nf'

format_ch = Channel.of( 'bed', 'tab')

workflow clean_cnv {
    take: 
    cnv
    pfb    
    genes
    links

    main:
    cnv
        | FILTER
        | combine(pfb)
        | CLEAN
        | groupTuple(by: 0)
        | combine(genes)
        | combine(links)
        | SCAN
        | combine(format_ch)
        | VISUALIZE

    emit:
    filtered    = FILTER.out
    cleaned     = CLEAN.out
    scanned     = SCAN.out
    visualized  = VISUALIZE.out
}

workflow {
    cnv = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    pfb = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    genes   = Channel.fromPath(params.refgene)
    links   = Channel.fromPath(params.reflink)
    
    clean_cnv(cnv, pfb, genes, links)
}