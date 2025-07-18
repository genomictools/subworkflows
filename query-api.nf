#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { PREDICT } from '../modules/predict.nf'
include { PLOT }    from '../modules/plot.nf'

workflow query_api {
    take:
    variants

    main:
    variants
        | PREDICT
        | filter { it.last().size() == 0 }
        | PLOT

    emit:
    predictions = PREDICT.out
    plots       = PLOT.out
}

// Workflow
workflow {
    // Define input from file
    variants_ch = Channel.fromPath(params.cohort_info)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.variant, row.ontology, row.assay, row.sequence_length ]}

    query_api(variants_ch)
}
