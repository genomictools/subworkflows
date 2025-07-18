#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { PREDICT } from '../modules/predict.nf'
include { PLOT }    from '../modules/plot.nf'

workflow predict_variant {
    take:
    variants

    main:
    // Predict variants using the Alphagenome API
    variants
        | PREDICT
        | set { predictions }

    // Generate plots if requested
    if ( params.plots ) {
        predictions
            | filter { it.last().size() == 0 }
            | PLOT
            | set { plots }
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

    predict_variant(variants_ch)
}
