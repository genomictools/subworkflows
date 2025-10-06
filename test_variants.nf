#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MATRIXEQTL }     from '../modules/matrixeqtl/matrixeqtl.nf'

tools_ch = Channel.of(params.tools.split(','))

workflow test_variants {
    take:
    cohorts
    tools

    main:
    cohorts | combine( tools ) | filter { it.last() == 'matrixeqtl' } | MATRIXEQTL
}

workflow  {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.category,
            file(row.snps), file(row.traits), file(row.covariates)
        ] }

    test_variants( cohorts_ch, tools_ch )
}
