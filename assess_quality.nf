#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ASSESS }    from '../modules/assess.nf'

workflow assess_quality {
    take: 
    cnv
    signal

    main:
    signal
        | combine(cnv, by: [0,1])
        | groupTuple(by: [0,1])
        | ASSESS

    emit:
    qc  = ASSESS.out
}

workflow {
    cnv    = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    signal = Channel.fromPath(params.signal) | map { [ it.simpleName, it ] }

    assess_quality(cnv, signal)
}