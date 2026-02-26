#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PLOT }        from '../modules/exomedepth/plot.nf'
include { VISUALIZE }   from '../modules/bamsignals/visualize.nf'

workflow visualize_cnv {
    take:
    calls
    controls
    coverage

    main:
    // Collect output in a file & filter
    calls
        | map { it[0,1,2,4]}
        | splitCsv(header: true, sep: '\t')
        | map { cohort, gene, key, cnv -> [ cohort: cohort, gene: gene, key: key ] + cnv }
        | filter { it['BF'].toDouble() >= params.BF }
        | filter { it['reads.ratio'] != 'NA' && it['reads.ratio'] != 'Inf'}
        | filter { it['reads.ratio'].toDouble() >= params.dup_ratio || it['reads.ratio'].toDouble() <= params.del_ratio}
        | map { [ it.cohort, it.gene, it.key ] }
        | set { filtered_calls }

    // Plot genes with CNVs
    calls
        | combine( filtered_calls, by: [0,1,2] )
        | PLOT
        | set { cnv_plots }

    // Visualize bams
    controls
        | splitCsv(header: true, sep: '\t')
        | map { cohort, key, ref -> [ cohort: cohort, key: key ] + ref }
        | filter { it.selected == 'TRUE' }
        | map { [it.cohort, [it['ref.samples'], it['test.samples']], it.key]}
        | transpose
        | unique
        | combine(coverage, by: [0,1])   // Combine by cohort and key
        | groupTuple(by: [0,2,4])              // Group by cohort, key and gene
        | map { it[0,4,2,1,6] }                // Extract cohort, gene, key, samples, and coverage
        | combine(filtered_calls, by: [0,1,2]) // Combine by cohort, gene and key
        | VISUALIZE
        | set { cnv_coverage }

    emit:
    cnv_plots
    cnv_coverage
}

workflow {
    calls_ch = Channel.fromPath(params.calls)
        | splitCsv(header: true, sep: '\t')
        | map { row -> [ row.cohort, row.gene, row.key, file(row.object), file(row.tsv), row.nmarker ] }

    controls_ch = Channel.fromPath(params.controls)
        | splitCsv(header: true, sep: '\t')
        | map { row -> [ row.cohort, row.key, file(row.ref) ] }

    coverage_ch = Channel.fromPath(params.coverage)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, row.gene, row.feature, file(row.depth) ] }

    visualize_cnv(calls_ch, controls_ch, coverage_ch )
}