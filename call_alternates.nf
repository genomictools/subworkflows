#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'
include { DETECT }    from '../modules/detect.nf'
include { EXCLUDE }   from '../modules/exclude.nf'
include { FILTER }    from '../modules/filter.nf'
include { CLEAN }     from '../modules/clean.nf'
include { ASSESS }    from '../modules/assess.nf'

type_ch     = Channel.of(params.type.split(','))
format_ch   = Channel.of(params.format.split(','))

workflow call_alternates {
    take: 
    gtc
    pfb
    gcm
    hmm
    hmm0
    exclude

    main:
    gtc
        | EXTRACT
        | filter { it.last().toInteger() > 1 }
        | ( params.adjust ? combine(gcm) : map { it } )
        | ( params.adjust ? ADJUST       : map { it } )
        | filter { it.last().toInteger() > 1 }
        | combine(pfb)
        | combine(hmm)
        | combine(hmm0)
        | combine(type_ch)
        | DETECT
        | filter { it.last().toInteger() > 1 }
        | filter { it[2] == 'cnv' }
        | ( params.exclude ? combine(exclude) : map { it } )
        | ( params.exclude ? EXCLUDE       : map { it } )
        | filter { it.last().toInteger() > 1 }
        | ( params.filter ? FILTER       : map { it } )
        | filter { it.last().toInteger() > 1 }
        | ( params.filter ? combine(pfb) : map { it } )
        | ( params.clean  ? CLEAN        : map { it } )
        | filter { it.last().toInteger() > 1 }
        | set { cnv }

    Channel.empty()
        | concat(DETECT.out)
        | filter { it[2] == 'loh' }
        | set { loh }

    EXTRACT.out
        | ( params.adjust ? concat(ADJUST.out) : map { it } )
        | set { signal }

    if ( params.assess ) {
        DETECT.out 
            | groupTuple(by: [0, 2])
            | ASSESS
    }

    emit:
    signal
    cnv
    loh
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }

    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm     = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }
    hmm     = Channel.fromPath(params.hmm)
    hmm0    = Channel.fromPath(params.hmm0)
    exclude = Channel.fromPath(params.exclude_regions)

    call_alternates(gtc_ch, pfb, gcm, hmm, hmm0, exclude)
}