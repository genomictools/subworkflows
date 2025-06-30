#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { VEP }         from '../modules/vep.nf'
include { CONCATINATE } from '../modules/concatinate.nf'

workflow run_vep {
    take:
    variants

    main:
    variants
        | VEP
        | groupTuple(by: [0,1,2])
        | CONCATINATE

    emit:
    CONCATINATE.out
}

// Workflow
workflow {
    // Define input from file
    variants_ch = Channel.fromPath(params.variants)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.id, file(row.file), file(row.index) ]}

    run_vep(variants_ch)
}
