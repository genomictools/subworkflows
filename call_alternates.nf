#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PENNCNV }   from '../modules/penncnv.nf'
include { QUANTISNP } from '../modules/quantisnp.nf'
include { RGADA }     from '../modules/rgada.nf'
include { CONVERT }   from '../modules/convert.nf'
include { COMBINE }   from '../modules/combine.nf'

workflow call_alternates {
    take: 
    signal
    pfb
    hmm
    levels
    type
    tools

    main:
    // PENNCNV
    tools
        | filter { it == 'penncnv' }
        | combine(type)
        | combine(hmm, by: 1)
        | map {[ it[1], it[0], it[2]]}
        | combine(pfb)
        | set { req }

    signal
        | filter { it[2] == 'adjusted' }
        | combine(req)
        | PENNCNV
        | filter { it.last().toInteger() > 1 }
        | set { penncnv }

    // QUANTISNP
    tools
        | filter { it == 'quantisnp' }
        | combine(levels, by: 0)
        | set { req }

    signal
        | filter { it[2] == 'merged' }
        | combine(req)
        | QUANTISNP
        | set { quantisnp }

    // RGADA
    tools
        | filter { it == 'rgada' }
        | set { req }

    signal
        | filter { it[2] == 'merged' }
        | combine(req)
        | RGADA
        | set { rgada }

    // Combine calls
    quantisnp
        | concat(rgada)
        | map {[
            it[0], it[1], it[2],
            it[3].split(',').toList(),
            it[4], it[5],
            it[6].split(',').toList()
        ]}
        | transpose
        | filter { it.last().toInteger() > 1 }
        | combine(pfb)
        | CONVERT
        | concat(penncnv)
        | groupTuple(by: [0, 2, 3])
        | COMBINE
        | set { calls }

    emit:
    calls
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }

    pfb     = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    levels  = Channel.fromPath(params.levels) | map { [ it.simpleName, it ] }
    type_ch = Channel.of(params.type.split(','))
    hmm     = Channel.empty()
        | ( params.hmm  != null ? concat(Channel.of(['cnv', file(params.hmm)]))  : Channel.empty() )
        | ( params.hmm0 != null ? concat(Channel.of(['loh', file(params.hmm0)])) : Channel.empty() )
        | combine(type_ch, by: 0)

    type_ch  = Channel.of(params.type.split(','))
    tools_ch = Channel.of(params.tools.split(','))

    call_alternates(gtc_ch, pfb, hmm, levels, type_ch, tools_ch)
}