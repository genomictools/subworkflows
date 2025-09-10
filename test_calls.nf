#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CCTEST }      from '../modules/cctest.nf'
include { FAMILY }      from '../modules/family.nf'
include { VALIDATE }    from '../modules/validate.nf'

workflow test_calls {
    take:
    signal 
    calls
    consensus
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
        | ( params.adjust ? filter { it[2] == 'adjusted' } : filter { it[2] == 'raw' } )
        | groupTuple(by: [0,2]) 
        | set { combined_signal}

    calls
        | combine(test)
        | filter { it.last() == 'family' }
        | map {[ it[0], it[2], it[1], it[3], it[6] ]}
        | combine(hmm, by: 1)
        | map {[ it[1], it[2], it[0], it[3], it[4], it[5] ]}
        | combine(pfb)
        | combine(pedigree, by: 0)
        | combine(combined_signal, by: 0)
        | FAMILY
        // | set { tested }
        | view

    // Validation
    consensus
        | combine(test)
        | filter { it.last() == 'validate' }
        | combine(hmm, by: 1)
        | map {[ it[1], it[0], it[3], it[6], it[7] ]}
        | combine(pfb)
        | combine(combined_signal, by: 0)
        | VALIDATE
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