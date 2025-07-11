#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT }     from '../modules/convert.nf'
include { PRUNE }       from '../modules/prune.nf'
include { COMBINE }     from '../modules/combine.nf'
include { TEST }        from '../modules/test.nf'
include { PLOT }        from '../modules/plot.nf'

test_ch      = Channel.of(params.tests.split(','))

workflow test_variants {
    take:
    genotypes
    
    main:
    genotypes
        | CONVERT
        | PRUNE
        | groupTuple(by: 0)
        | COMBINE
        | combine(test_ch)
        | TEST
        | PLOT

    emit:
    tests = TEST.out
    plots = PLOT.out
}

workflow  {
    genotypes_ch = Channel.fromPath(params.genotypes)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.key,
            file(row.file), file(row.index), file(row.samples)
        ] }

    test_variants(genotypes_ch)
}
