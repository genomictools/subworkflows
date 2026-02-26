#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SELECT }      from '../modules/exomedepth/select.nf'
include { CALL }        from '../modules/exomedepth/call.nf'
include { REPORT }      from '../modules/rocker/report.nf'

workflow call_exome_depth {
    take:
    counts
    coverage

    main:
    counts
        | branch {
            ref_controls : it[2] == 'control' & it[3] == 'reference'
            ref_cases    : it[2] != 'control' & it[3] == 'reference'
            gen_controls : it[2] == 'control' & it[3] != 'reference'
            gen_cases    : it[2] != 'control' & it[3] != 'reference'
        }
        | set { counts }

    counts.ref_cases
        | combine(counts.ref_controls, by: [0,3,4]) // Combine by cohort, range, feature
        | groupTuple(by: [0,3,5])                   // Group by cohort, key and counts
        | map { it[0,3,5,6,8] }                     // Extract cohort, key, counts, ref, and ref counts
        | SELECT
        | set { controls }

    controls
        | collectFile(
            storeDir: "${params.output_dir}/select",
            keepHeader: true,
            skip: 1
        ) { [ "${it[0]}.select.tsv", it[2]] }
        | map { file -> [ file.name.split('\\.')[0], 'samples', file, file.readLines().size()] }
        | filter { it.last() > 1 }
        | REPORT
        | set { reports }

    // Aggregate reference ranges, in each case sample
    controls
        | splitCsv(header: true, sep: '\t')
        | map { cohort, key, ref -> [ cohort: cohort, key: key ] + ref }
        | filter { it.selected == 'TRUE' }
        | map { [ it.cohort, it.key, it['ref.samples']] }
        | combine( counts.gen_cases, by: [0,1] )
        | map { cohort, key, ref_key, type, gene, feature, count ->
            [ cohort, ref_key, type, gene, count, key ]
        }
        | combine( counts.gen_controls, by: [0,1,3] ) // Combine by cohort, ref_key and gene
        | groupTuple(by: [0,2,5,4])              // Group by cohort, gene, key and counts
        | map { it[0,2,5,4,8] }                  // Extract cohort, gene, key, counts, and ref counts
        | CALL
        | filter { it.last().toInteger() > 1 }
        | set { calls }

    // Save output to file
    calls
        | collectFile(
            storeDir: "${params.output_dir}/calls",
            keepHeader: true,
            skip: 1
        ) { [ "${it[0]}.cnv.tsv", it[4]] }
        | map { file -> [ file.name.split('\\.')[0], 'calls', file, file.readLines().size()] }
        | filter { it.last() > 1 }
        // | REPORT
        | set { reports }

    emit:
    calls    = calls
    controls = controls
    // reports = reports
}

workflow {
    counts_ch = Channel.fromPath(params.counts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, row.gene, row.feature, file(row.depth) ] }

    coverage_ch = Channel.fromPath(params.coverage)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, row.gene, row.feature, file(row.depth) ] }

    call_exome_depth(counts_ch, coverage_ch )
}