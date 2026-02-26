#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CALL }   from '../modules/coveragemaster/call.nf'
include { SELECT } from '../modules/coveragemaster/select.nf'

workflow call_coverage_master {
    take:
    coverage
    stats

    main:
    coverage
        | filter { it[3] != 'reference' }
        | combine(stats, by: [0,1,2])
        | branch {
            controls : it[2] == 'control'
            cases    : it[2] != 'control'
        }
        | set { coverage }

    coverage.cases
        | combine(coverage.controls, by: [0,3])
        | groupTuple(by: [0,1,2,5,7]) 
        | map { it[0,1,2,5,7,11,13] }
        | SELECT
        | CALL
        | set { calls }

    emit:
    calls
}

workflow {
    coverage_ch = Channel.fromPath(params.coverage)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, row.gene, row.feature, file(row.depth) ] }

    stats_ch = Channel.fromPath(params.stats)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, file(row.stats) ] }

    call_coverage_master( coverage_ch, stats_ch)
}