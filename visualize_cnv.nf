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
    cnv
    genes
    links
    signal
    pfb

    main:
    cnv
        | groupTuple(by: [0, 2])
        | combine(genes)
        | combine(links)
        | combine(features_ch)
        | SCAN
        | filter { it.last().toInteger() > 1 }
        | set { annotated }
    
    if ( params.export ) {
    annotated
        | combine(format_ch)
        | EXPORT
        | set { tables }
    } else {
        Channel.empty() | set { tables }
    }

    if ( params.plot ) {
    signal
        | ( params.adjust ? filter { it[2] == 'adjusted' } : map { it } )
        | groupTuple(by: [ 0, 2 ])
        | set { plot_signal }

    annotated
        | filter { it[1] == 'gene' }
        | combine(plot_signal, by: 0)
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
    cnv     = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    genes   = Channel.fromPath(params.refgene)
    links   = Channel.fromPath(params.reflink)

    signal  = Channel.fromPath(params.signal)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }
    
    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }

    visualize_cnv(cnv, genes, links, signal, pfb)
}