#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { ALPHAGENOME } from '../modules/alphagenome.nf'
include { RESHAPE }     from '../modules/reshape.nf'
include { FORMAT }      from '../modules/format.nf'
include { CONCATINATE } from '../modules/concatinate.nf'

workflow run_alphagenome {
    take:
    variants

    main:
    variants
        | ALPHAGENOME
        | FORMAT
        | RESHAPE
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

    run_alphagenome(variants_ch)
}
