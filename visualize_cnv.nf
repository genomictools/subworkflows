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
        | set { annotations }

    // export as tables
    format_ch = Channel.of(params.format.split(','))
    Channel.empty() | set { tables }
    if ( params.export ) {
    annotations
        | combine(format_ch)
        | EXPORT
        | set { tables }
    }

    // heatmpas
    Channel.empty() | set { heatmaps }
    if ( params.heatmap ) {
    annotations
        | filter { it.last().toInteger() > 1 }
        | HEATMAP
        | set { heatmaps }
    }

    // scatter plots
    Channel.empty() | set { scatter_plots }
    if ( params.scatter ) {
    // get samples and features
    annotations
        | splitCsv(elem: 4, header: false, strip: true, sep: "\t")
        | map { cohort, tool, feature, type, row, log, nmarkers -> 
            def cnv    = row[0].replaceAll(' +', '\t').split('\t')
            def region = cnv[0].replaceAll(':|-', '\\_')
            def sample = cnv[4]

            if      ( feature == 'refgene' ) { select_list = params.genelist } 
            else if ( feature == 'anno' )    { select_list = params.bandlist }
            def select_list  = file(select_list).readLines()

            def feature_list = row[1].split(',').toList().intersect(select_list)
            def n_features   = feature_list.size()

            return [ cohort, sample, type, tool, region, feature, feature_list, n_features ]
        }
        | filter { it.last().toInteger() > 0 }
        | combine(signal, by: [0,1])
        | ( params.adjust ? filter { it[8] == 'adjusted' } : filter { it[8] == 'raw' } )
        | combine(pfb)
        | SCATTER
        | set { scatters }
    }

    emit:
    annotations
    tables
    heatmaps
    scatters
}

workflow {
    calls   = Channel.fromPath(params.calls) | map { [ it.simpleName, it ] }
    signal  = Channel.fromPath(params.signal)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }
    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }

    visualize_cnv(calls, pfb, signal)
}
