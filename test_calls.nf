#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CCTEST }      from '../modules/penncnv/cctest.nf'
include { FAMILY }      from '../modules/penncnv/family.nf'
include { VALIDATE }    from '../modules/penncnv/validate.nf'

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
    // Return cohort, family name, family size and test
    pedigree
        | map { 
            def lines = it[1].toFile().text.split('\n').findAll { !it.trim().isEmpty() }
            def familyGroups = lines.groupBy { it.split(/\s+/)[0] }
            
            familyGroups.collect { fam, members ->
                def size = members.size()
                def type = fam == '0' ? 'cctest' : 
                          size < 3 ? 'cctest' : 
                          size == 3 ? 'trio' : 
                          size == 4 ? 'quartet' : 'family'
                [it[0], fam, size, type]
            }
        }
        | flatMap()
        | groupTuple(by: [0,3])
        | branch {
            cc  : it.last() == 'cctest'
            fam : it.last() != 'cctest'
        }
        | set { test }

    calls
        | combine(test.cc, by: 0)
        | combine(pedigree, by: 0)
        | combine(pfb)
        | CCTEST
        | set { cc_test }

    // Combined signal
    signal
        | ( params.adjust ? filter { it[2] == 'adjusted' } : filter { it[2] == 'raw' } )
        | groupTuple(by: [0,2]) 
        | set { combined_signal }

    calls
        | combine(test.fam, by: 0)
        | combine(pedigree, by: 0)
        | combine(pfb)
        | map { item -> [0, 2, 1, *(3..<item.size())].collect { item[it] } }
        | combine(hmm, by: 1)
        | map { item -> [1, 2, 0, *(3..<item.size())].collect { item[it] } }
        | combine(combined_signal, by: 0)
        | FAMILY
        | set { family_test }

    // Validation
    consensus
        | combine(Channel.of ("validate"))
        | combine(pfb)
        | combine(hmm, by: 1)
        | map { item -> [1, 0, *(2..<item.size())].collect { item[it] } }
        | combine(combined_signal, by: 0)
        | VALIDATE
        | set { validation_test }

    // emit:
    // tested = tested
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