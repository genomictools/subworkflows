#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ANNOTATE }  from '../modules/penncnv/annotate.nf'
include { EXPORT }    from '../modules/penncnv/export.nf'
include { HEATMAP }   from '../modules/cnvr/heatmap.nf'
include { SCATTER }   from '../modules/cnvr/scatter.nf'

workflow visualize_cnv {
    take: 
    calls
    pfb
    signal

    main:
    // annotate calls
    features = Channel.of (
            [ 'refgene', params.refgene ], 
            [ 'refexon', params.refexon ], 
            [ 'anno', params.anno ] 
        )
        | filter { it.last() != null }
        | map { [it.first(), file(it.last())] }

    calls
        | combine(features)
        | ANNOTATE
        | filter { it.last().toInteger() > 1 }
        | branch { 
            gene     : it[2] == 'refgene'
            segments : it[2] == 'anno'
         }
        | set { annotated }
    
    // export as tables
    format_ch = Channel.of(params.format.split(','))
    Channel.empty() | set { tables }
    if ( params.export ) {
    annotated.gene
        | combine(format_ch)
        | EXPORT
        | set { tables }
    }

    // plots
    Channel.fromPath(params.genelist)
        | map { ['refgene', it.readLines()] }
        | set { genelist_ch }

    genelist_ch
        | transpose
        | set { genelist_ch_t }

    // heatmpas
    Channel.empty() | set { heatmaps }
    if ( params.heatmap ) {
    annotated.gene
        | filter { it.last().toInteger() > 1 }
        | combine(genelist_ch)
        | HEATMAP
        | set { heatmaps }
    }

    // scatter plots
    Channel.empty() | set { scatter_plots }
    if ( params.scatter ) {
    // get samples and features
    annotated.gene
        | splitCsv(elem: 4, header: false, strip: true, sep: "\t")
        | multiMap { cohort, tool, feature, type, row, log, nmarkers -> 
            def cnv    = row[0].replaceAll(' +', '\t').split('\t')
            def region = cnv[0].replaceAll(':|-', '\\_')
            def sample = cnv[4].split('\\.')[1]
            def gene   = row[1].split(',').toList()
            samples  : [ cohort, sample, type, tool, region ]
            features : [ feature, gene, region ]
        }
        | set { annotations }

    // subset feagures
    annotations.features
        | transpose
        | combine(genelist_ch_t, by: [0,1])
        | unique
        | groupTuple(by: [0,2])
        | map { tuple(it[-1], *it[0..it.size()-2]) }
        | set { anno_features }

    // combine signal with annotations, and plot
    signal
        | ( params.adjust ? filter { it[2] == 'adjusted' } : filter { it[2] == 'raw' } )
        | combine(annotations.samples, by: [0,1])
        | map { tuple(it[-1], *it[0..it.size()-2]) }
        | combine(anno_features, by: 0)
        | map { tuple(*it[1..it.size()-1], it[0]) }
        | combine(pfb)
        | SCATTER
        | set { scatter_plots }
    }

    emit:
    segments = annotated.segments
    genes    = annotated.gene
    tables
    heatmaps
    scatter_plots
}

workflow {
    calls   = Channel.fromPath(params.calls) | map { [ it.simpleName, it ] }
    signal  = Channel.fromPath(params.signal)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }
    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }

    visualize_cnv(calls, pfb, signal)
}
