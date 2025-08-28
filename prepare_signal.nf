#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'
include { MERGE }     from '../modules/merge.nf'

workflow prepare_signal {
    take: 
    gtc
    pfb
    gcm

    main:
    // Extract and adjust signal
    gtc
        | EXTRACT
        | filter { it.last().toInteger() > 1 }
        | ( params.adjust ? combine(gcm) : map { it } )
        | ( params.adjust ? ADJUST       : map { it } )
        | filter { it.last().toInteger() > 1 }
        | ( params.merge ? combine(pfb) : map { it } )
        | ( params.merge ? MERGE        : map { it } )

    EXTRACT.out
        | ( params.adjust ? concat(ADJUST.out) : map { it } )
        | ( params.adjust ? concat(MERGE.out)  : map { it } )
        | set { signal }

    emit:
    signal
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',§')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }
    pfb    = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm    = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }

    prepare_signal(gtc_ch, pfb, gcm)
}