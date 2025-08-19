#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SCAN }      from '../modules/scan.nf'
include { EXPORT }    from '../modules/export.nf'
include { PLOT }      from '../modules/plot.nf'

format_ch   = Channel.of(params.format.split(','))
features_ch = Channel.of(params.features.split(','))
plot_type_ch= Channel.of(params.plot_type.split(','))

workflow visualize_cnv {
    take: 
    calls
    pfb
    signal
    genes
    links


    main:
    // annotate calls
    calls
        | combine(genes)
        | combine(links)
        | combine(features_ch)
        | SCAN
        | filter { it.last().toInteger() > 1 }
        | set { annotated }
    
    // export as tables
    if ( params.export ) {
    annotated
        | combine(format_ch)
        | EXPORT
        | set { tables }
    } else {
        Channel.empty() | set { tables }
    }

    // plot calls
    if ( params.plot ) {
    annotated
        | combine(signal, by: 0)
        | groupTuple(by: [ 0,1,2,3,4,7 ])
        | combine(pfb)
        | combine(plot_type_ch)
        | PLOT
        | set { plots }
    } else {
        Channel.empty() | set { plots }
    }

    emit:
    annotated
    tables
    plots
}

workflow {
    calls   = Channel.fromPath(params.calls) | map { [ it.simpleName, it ] }
    signal  = Channel.fromPath(params.signal)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }
    
    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    genes   = Channel.fromPath(params.refgene)
    links   = Channel.fromPath(params.reflink)

    visualize_cnv(calls, pfb, signal, genes, links)
}