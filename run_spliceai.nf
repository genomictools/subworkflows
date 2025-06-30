#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { SPLICEAI }    from '../modules/spliceai.nf'
include { CONCATINATE } from '../modules/concatinate.nf'

workflow run_spliceai {
    take:
    variants

    main:
    variants
        | SPLICEAI
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

    run_spliceai(variants_ch)
}
