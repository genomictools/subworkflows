#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TABULATE }     from '../modules/bcftools/tabulate.nf'
include { SUMMARIZE }    from '../modules/rocker/summarize.nf'

workflow summarize_controls {
    take:
    genotypes
    annotations

    main:
    genotypes
        | TABULATE
        | groupTuple(by: [0,2])
        | SUMMARIZE
        | set { summary }

    emit:
    summary = summary
}

workflow  {
    genotypes_ch = Channel.fromPath(params.genotypes)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.category, row.file, row.index, row.n_vars ] }

    annotations_ch = Channel.fromPath(params.annotations)
        | map { row -> [ row.cohort, row.key, row.category, row.variable, row.file ] }

    summarize_controls( genotypes_ch, annotations_ch )
}
