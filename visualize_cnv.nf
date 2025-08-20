#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ANNOTATE }  from '../modules/annotate.nf'
include { EXPORT }    from '../modules/export.nf'
include { HEATMAP }   from '../modules/heatmap.nf'
include { SCATTER }   from '../modules/scatter.nf'

format_ch   = Channel.of(params.format.split(','))
features_ch = Channel.of(params.features.split(','))

workflow visualize_cnv {
    take: 
    calls
    pfb
    signal
    genes
    links
    genelist

    main:
    // annotate calls
    calls
        | combine(genes)
        | combine(links)
        | combine(features_ch)
        | ANNOTATE
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
    cohort_gene = genelist | groupTuple(by: 0)
    if ( params.heatmap ) {
    annotated
        | combine(cohort_gene, by: 0)
        | HEATMAP
        | set { plots }
    } else {
        Channel.empty() | set { plots }
    }

    if ( params.scatter ) {
    annotated
        | splitCsv(elem: 3, header: false, strip: true, sep: "\t")
        | map { cohort, feature, type, row, log, nmarkers ->
            def gene = row[1].split(',').toList()
            return [cohort, gene, row[0] ]
        }
        | transpose
        | unique()
        | groupTuple(by: [0,1])
        | ( params.genelist != null ? combine(genelist, by: [0,1]) : map { it })
        | combine(signal, by: 0)
        | combine(pfb)
        | SCATTER
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

    genelist_ch = Channel.empty()
        | ( params.genelist != null ? concat(Channel.of(file(params.genelist))) : Channel.empty() )
        | splitCsv(header: true)
        | map { row -> [ row.cohort, row.gene ] }

    visualize_cnv(calls, pfb, signal, genes, links, genelist_ch)
}
