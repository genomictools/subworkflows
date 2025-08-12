#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SCAN }      from '../modules/scan.nf'
include { VISUALIZE } from '../modules/visualize.nf'

format_ch   = Channel.of('bed', 'tab')
features_ch = Channel.of(params.features.split(','))

workflow visualize_cnv {
    take: 
    cnv
    genes
    links

    main:
    cnv
        | groupTuple(by: [0, 2])
        | combine(genes)
        | combine(links)
        | combine(features_ch)
        | SCAN
        | combine(format_ch)
        | VISUALIZE

    emit:
    scanned     = SCAN.out
    visualized  = VISUALIZE.out
}

workflow {
    cnv = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    genes   = Channel.fromPath(params.refgene)
    links   = Channel.fromPath(params.reflink)
    
    visualize_cnv(cnv, genes, links)
}