#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXCLUDE }   from '../modules/exclude.nf'
include { FILTER }    from '../modules/filter.nf'
include { CLEAN }     from '../modules/clean.nf'
include { ASSESS }    from '../modules/assess.nf'

workflow clean_calls {
    take: 
    cnv
    pfb
    exclude

    main:
    cnv
        | ( params.exclude ? combine(exclude) : map { it } )
        | ( params.exclude ? EXCLUDE          : map { it } )
        | filter { it.last().toInteger() > 1 }
        | ( params.filter ? FILTER       : map { it } )
        | filter { it.last().toInteger() > 1 }
        | ( params.filter ? combine(pfb) : map { it } )
        | ( params.clean  ? CLEAN        : map { it } )
        | filter { it.last().toInteger() > 1 }
        | set { cleaned }
    
    // if ( params.assess ) {
    //     DETECT.out 
    //         | groupTuple(by: [0, 2])
    //         | ASSESS
    //         | set { assessed }
    // }

    emit:
    cnv = cleaned
    // assessed
}

workflow {
    cnv     = Channel.fromPath(params.cnv) | map { [ it.simpleName, it ] }
    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    exclude = Channel.fromPath(params.exclude_regions)

    clean_calls(cnv, pfb, exclude)
}