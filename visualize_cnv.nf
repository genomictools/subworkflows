#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SCAN }      from '../modules/scan.nf'
include { VISUALIZE } from '../modules/visualize.nf'

format_ch   = Channel.of(params.format.split(','))
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
        | set { annotated }
    
    annotated
        | filter { it.last().toInteger() > 1 }
        | combine(format_ch)
        | VISUALIZE
        | filter { it.last().toInteger() > 1 }
        | set { tables }

    emit:
    annotated
    tables
}

workflow {
    cnv     = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    genes   = Channel.fromPath(params.refgene)
    links   = Channel.fromPath(params.reflink)
    
    visualize_cnv(cnv, genes, links)
}