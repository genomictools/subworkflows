#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT }     from '../modules/convert.nf'
include { PRUNE }       from '../modules/prune.nf'
include { COMBINE }     from '../modules/combine.nf'
include { TEST }        from '../modules/test.nf'
include { PLOT }        from '../modules/plot.nf'

test_ch = Channel.of(params.tests.split(','))

workflow test_variants {
    take:
    genotypes

    main:
    genotypes
        | CONVERT
        | ( params.prune ? PRUNE : map { it } )
        | groupTuple(by: [0, 2])
        | COMBINE
        | combine(test_ch)
        | TEST
        | ( params.plot ? PLOT : map { it })

    emit:
    tests = TEST.out
}

workflow  {
    genotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort,row.key,row.category,file(row.file),file(row.index),
            row.n_samples,row.n_variants,
            file(row.phenotype)
        ] }

    test_variants( genotypes_ch )
}
