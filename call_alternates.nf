#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT }   from '../modules/penncnv/convert.nf'
include { REPORT }    from '../modules/penncnv/report.nf'
include { PENNCNV }   from '../modules/penncnv/penncnv.nf'
include { QUANTISNP } from '../modules/quantisnp/quantisnp.nf'
include { RGADA }     from '../modules/rgada/rgada.nf'
include { SPLIT }     from '../modules/rocker/split.nf'

workflow call_alternates {
    take: 
    signal
    pfb

    main:
    tools = Channel.of(params.tools.split(','))
    // PENNCNV
    signal
        | ( params.adjust ? filter { it[2] == 'adjusted' } : filter { it[2] == 'raw' } )
        | combine(tools) | filter { it.last() == 'penncnv' }
        | combine(pfb)
        | PENNCNV
        | filter { it.last().toInteger() > 0 }
        | multiMap { cohort, key, tool, cnv, log, nmarker ->
            cnv  : [ cohort, key, tool, cnv ]
            log  : [ cohort, key, tool, log ]
            size : [ cohort, key, tool, nmarker ]
        }
        | set { penncnv }
    penncnv.cnv 
        | collectFile(storeDir: "${params.output_dir}/penncnv") { [ "${it[0]}.${it[2]}.combined.cnv", it[3]] }
        | map { file -> [file.name.split('\\.')[0], file.name.split('\\.')[1], file ] }
        | set { penncnv_cnv }
    penncnv.log 
        | collectFile(storeDir: "${params.output_dir}/penncnv") { [ "${it[0]}.${it[2]}.combined.log", it[3]] }
        | map { file -> [file.name.split('\\.')[0], file.name.split('\\.')[1], file ] }
        | set { penncnv_log }
    penncnv.size 
        | groupTuple(by: [0, 2])
        | map { cohort, key, tool, nmarker -> 
            def sum = nmarker.collect { it.toInteger() }.sum()
            [cohort, tool, sum]
        }
        | set { penncnv_size }

    penncnv_cnv
        | combine(penncnv_log,  by: [0, 1])
        | combine(penncnv_size, by: [0, 1])
        | unique
        | set { penncnv }

    // QUANTISNP
    signal
        | filter { it[2] == 'merged' }
        | combine(tools) | filter { it.last() == 'quantisnp' }
        | QUANTISNP
        | set { quantisnp }

    // RGADA
    signal
        | filter { it[2] == 'merged' }
        | combine(tools) | filter { it.last() == 'rgada' }
        | RGADA
        | set { rgada }

    // Combine calls
    quantisnp
        | concat(rgada)
        | collectFile(keepHeader: true, storeDir: "${params.output_dir}/") { [ "${it[2]}/${it[0]}.${it[2]}.combined.cnv", it[3]] }
        | map { file -> [file.name.split('\\.')[0], file.name.split('\\.')[1], file ] }
        | combine(pfb)
        | CONVERT
        | concat(penncnv)
        | SPLIT
        | map { cohort, tool, type, cnv, log, nmarker -> [
            cohort, tool,
            type.split(',').toList(),
            cnv, log,
            nmarker.split(',').toList()
        ]}
        | transpose
        | filter { it.last().toInteger() > 1 }
        | set { calls }

    // QC Reports
    reports = Channel.empty()
    if ( params.report ) {
    penncnv
        | ( params.report ? REPORT : map { it } )
        | set { reports }
    }

    emit:
    calls
    reports
}

workflow {
    signal_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }

    pfb_ch      = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }

    call_alternates(signal_ch, pfb_ch)
}
