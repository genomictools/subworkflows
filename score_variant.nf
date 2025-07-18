#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { SCORE }   from '../modules/score.nf'

workflow score_variant {
    take:
    variants

    main:
    // Download recommended scores if requested
    variants
        | SCORE
        | set { scores }

    emit:
    scores
}

// Workflow
workflow {
    // Define input from file
    variants_ch = Channel.fromPath(params.cohort_info)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.variant, row.organism, row.sequence_length ]}
        | unique

    score_variant(variants_ch)
}
