#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'
include { DETECT }    from '../modules/detect.nf'

workflow call_alternates {
    take: 
    gtc
    pfb
    gcm
    hmm

    main:
    // Extract and adjust signal
    gtc
        | EXTRACT
        | filter { it.last().toInteger() > 1 }
        | ( params.adjust ? combine(gcm) : map { it } )
        | ( params.adjust ? ADJUST       : map { it } )
        | filter { it.last().toInteger() > 1 }
        | set { signal }

    // call alternates
    signal
        | combine(hmm)
        | combine(pfb)
        | DETECT
        | filter { it.last().toInteger() > 1 }
        | set { calls }

    emit:
    signal
    calls
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }

    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm     = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }
    type_ch = Channel.of(params.type.split(','))
    hmm     = Channel.empty()
        | ( params.hmm  != null ? concat(Channel.of(['cnv', file(params.hmm)]))  : Channel.empty() )
        | ( params.hmm0 != null ? concat(Channel.of(['loh', file(params.hmm0)])) : Channel.empty() )
        | combine(type_ch, by: 0)

    call_alternates(gtc_ch, pfb, gcm, hmm)
}