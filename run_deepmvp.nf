#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { CSQ }         from '../modules/csq.nf'
include { DEEPMVP }     from '../modules/deepmvp.nf'
include { FORMAT }      from '../modules/format.nf'
include { RESHAPE }     from '../modules/reshape.nf'
include { CONCATINATE } from '../modules/concatinate.nf'

workflow run_deepmvp {
    take:
    variants

    main:
    variants
        | CSQ
        | filter { it.last().toInteger() > 1 }
        | DEEPMVP
        | filter { it.last().toInteger() > 1 }
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

    run_deepmvp(variants_ch)
}
