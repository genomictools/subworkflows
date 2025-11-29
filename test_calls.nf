#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ROH }         from '../modules/plink/roh.nf'
include { CCTEST }      from '../modules/penncnv/cctest.nf'
include { FAMILY }      from '../modules/penncnv/family.nf'
include { VALIDATE }    from '../modules/penncnv/validate.nf'

workflow test_calls {
    take:
    signal
    genotypes
    calls
    consensus
    pedigree
    pfb

    main:
    // ROH
    genotypes
        | ROH
        | set { roh }

    // TEST
    pedigree
        | flatMap { cohort, file ->
            def lines = file.text.readLines().findAll { it.trim() }
            def familyGroups = lines.groupBy { it.split(/\s+/)[0] }
            familyGroups.collect { fam, members ->
                def sample_ids = members.collect { it.split(/\s+/)[1] }
                def size = sample_ids.size()
                def type = (fam == '0' || size < 3) ? 'cctest'
                          : size == 3 ? 'trio'
                          : size == 4 ? 'quartet'
                          : 'family'
                def famid = size >= 3 ? fam : '0'
                [ cohort, famid, type, sample_ids ]
            }
        }
        | groupTuple(by: [0,1,2])
        | map { cohort, famid, type, sample_ids -> [ cohort, sample_ids.flatten(), sample_ids.flatten().size(), famid, type ] }
        | branch {
            cc  : it.last() == 'cctest'
            fam : it.last() != 'cctest'
        }
        | set { test }

    // CCTEST
    calls
        | combine(test.cc, by: 0)
        | combine(pedigree, by: 0)
        | combine(pfb)
        | CCTEST
        | set { cc_test }

    // FAMILY
    test.fam
        | transpose
        | combine(signal, by: [0,1])
        | ( params.adjust ? filter { it[5] == 'adjusted' } : filter { it[5] == 'raw' } )
        | groupTuple(by: [0,2,3,4,5])
        | set { family_signal }

    calls
        | combine(family_signal, by: 0)
        | combine(pedigree, by: 0)
        | combine(pfb)
        | FAMILY
        | set { family_test }

    // VALIDATE
    signal
        | ( params.adjust ? filter { it[2] == 'adjusted' } : filter { it[2] == 'raw' } )
        | groupTuple(by: [0,2])
        | set { combined_signal }

    consensus
        | combine(Channel.of ("validate"))
        | combine(combined_signal, by: 0)
        | combine(pfb)
        | VALIDATE
        | set { validation_test }

    Channel.empty() 
        | concat(cc_test)
        | concat(family_test)
        | concat(validation_test)
        | set { tested }

    emit:
    tested = tested
    roh    = roh
}

workflow {
    signal  = Channel.fromPath(params.signal)   | map { [ it.simpleName, it ] }
    genotypes  = Channel.fromPath(params.genotypes)   | map { [ it.simpleName, it ] }
    calls   = Channel.fromPath(params.cnv)      | map { [ it.simpleName, it ] }
    pedigree= Channel.fromPath(params.pedigree) | map { [ it.simpleName, it ] }
    pfb     = Channel.fromPath(params.pfb)      | map { [ it.simpleName, it ] }

    test_calls(signal, genotypes, calls, pedigree, pfb)
}