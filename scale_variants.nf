#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXCLUDE }     from '../modules/exclude.nf'
include { PRUNE }       from '../modules/prune.nf'
include { SCALE }       from '../modules/scale.nf'

workflow scale_variants {
    take:
    genotypes

    main:
    genotypes
        | ( params.prune  ? PRUNE  : map { it } )
        | filter { it.last().toInteger() > 0 }
        | ( params.exclude ? EXCLUDE : map { it } )
        | filter { it.last().toInteger() > 0 }
        | SCALE
        | map { [it[0], it[1], it[2]]}
        | set { scaled }

    emit:
    scaled
}

workflow  {
    genotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.category,
            file(row.file), file(row.index),
            row.n_samples,row.n_variants
        ] }

    scale_variants( genotypes_ch )
}