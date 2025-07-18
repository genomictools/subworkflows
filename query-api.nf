#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { PREDICT } from '../modules/predict.nf'
include { PLOT }    from '../modules/plot.nf'
include { SCORE }   from '../modules/score.nf'

workflow query_api {
    take:
    variants

    main:
    // Predict variants using the Alphagenome API
    variants
        | PREDICT
        | set { predictions }

    // Generate plots if requested
    Channel.empty() | set { scores }
    if ( params.plots ) {
        predictions
            | filter { it.last().size() == 0 }
            | PLOT
            | set { plots }
    }

    // Download recommended scores if requested
    Channel.empty() | set { scores }
    if ( params.scores ) {
        variants
            | map { it -> [ it[0], it[1], it[2], it[5] ] }
            | unique
            | SCORE
            | set { scores }
    }

    emit:
    predictions
}

// Workflow
workflow {
    // Define input from file
    variants_ch = Channel.fromPath(params.cohort_info)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.variant, row.organism, row.ontology, row.assay, row.sequence_length ]}

    query_api(variants_ch)
}
