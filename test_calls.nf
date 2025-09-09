#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CCTEST }      from '../modules/cctest.nf'
include { FAMILY }      from '../modules/family.nf'

workflow test_calls {
    take:
    signal 
    calls
    pedigree
    pfb
    hmm 
    test

    main:
    // CCTEST
    calls
        | combine(test)
        | filter { it.last() == 'cctest' }
        | combine(pedigree, by: 0)
        | combine(pfb)
        | CCTEST
        | set { tested }

    // FAMILY
    signal
        | groupTuple(by: [0,2]) 
        | set { combined_signal}

    calls
        | map {[ it[0], it[2], it[1] ]}
        | combine(hmm, by: 1)
        | map {[ it[1], it[2], it[0], it[3] ]}
        | combine(pfb)
        | combine(pedigree, by: 0)
        | combine(combined_signal, by: 0)
        | set { req }

    calls
        | combine(test)
        | filter { it.last() == 'family' }
        | combine(req, by: [0,1,2])
        | FAMILY
        | set { tested }

    emit:
    tested = tested
}

workflow {
    signal  = Channel.fromPath(params.signal)   | map { [ it.simpleName, it ] }
    calls   = Channel.fromPath(params.cnv)      | map { [ it.simpleName, it ] }
    pedigree= Channel.fromPath(params.pedigree) | map { [ it.simpleName, it ] }
    pfb     = Channel.fromPath(params.pfb)      | map { [ it.simpleName, it ] }
    hmm     = Channel.fromPath(params.hmm)      | map { [ it.simpleName, it ] }
    tests   = Channel.of(params.tests.split(','))

    test_calls(signal, calls, pedigree, pfb, hmm, tests)
}