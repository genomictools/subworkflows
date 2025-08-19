#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'
include { DETECT }    from '../modules/detect.nf'

type_ch     = Channel.of(params.type.split(','))

workflow call_alternates {
    take: 
    gtc
    pfb
    gcm
    hmm
    hmm0

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
        | combine(pfb)
        | combine(hmm)
        | combine(hmm0)
        | combine(type_ch)
        | DETECT
        | filter { it.last().toInteger() > 1 }
        | branch { 
            cnv : it[2] == 'cnv'
            loh : it[2] == 'loh'
        }
        | set { calls }

    emit:
    signal
    cnv = calls.cnv
    loh = calls.loh
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }

    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm     = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }
    hmm     = Channel.fromPath(params.hmm)
    hmm0    = Channel.fromPath(params.hmm0)

    call_alternates(gtc_ch, pfb, gcm, hmm, hmm0)
}