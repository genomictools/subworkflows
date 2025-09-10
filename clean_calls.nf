#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FILTER }    from '../modules/filter.nf'
include { EXCLUDE }   from '../modules/exclude.nf'
include { CLEAN }     from '../modules/clean.nf'
include { CNVR }      from '../modules/cnvr.nf'

workflow clean_calls {
    take: 
    calls
    pfb
    exclude
    cohort_size

    main:
    calls
        | filter { it[2] == 'cnv' } // TODO: remove when other types are supported
        | ( params.filter ? FILTER       : map { it } )
        | filter { it.last().toInteger() > 0 }
        | ( params.exclude ? combine(exclude) : map { it } )
        | ( params.exclude ? EXCLUDE          : map { it } )
        | filter { it.last().toInteger() > 0 }
        | ( params.clean ? combine(pfb) : map { it } )
        | ( params.clean  ? CLEAN        : map { it } )
        | filter { it.last().toInteger() > 0 }
        | set { cleaned }
    
    consensus = Channel.empty()
    if (params.consensus) {
    cleaned
        | groupTuple(by: [0,2])
        | combine(cohort_size, by: 0)
        | CNVR
        | set { consensus }
    }

    emit:
    calls = cleaned
    consensus = consensus
}

workflow {
    calls   = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    exclude = Channel.fromPath(params.exclude_regions)
    cohort_size = calls.count().map { it as Integer }
    clean_calls(calls, pfb, exclude, cohort_size)
}