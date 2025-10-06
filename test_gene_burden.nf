#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { READFREQUENCY }   from '../modules/rocker/readfrequency.nf'
include { BURDENTEST }      from '../modules/rocker/burdentest.nf'

model_ch = Channel.of('DOM', 'REC')

workflow test_gene_burden {
    take:
    cohorts
    model

    main:
    // Load cohorts
    cohorts
        | READFREQUENCY
        | branch {
            cases: it[1] == 'cases'
            controls: it[1] == 'controls'
        }
        | set { counts }

    // Run tests
    counts.cases 
        | combine(counts.controls, by: 2)
        | combine(model)
        | filter { it[0] != 'ALL' }
        | BURDENTEST
        | collectFile (
            keepHeader: true,
            storeDir: "${params.output_dir}/summary",
        )
        { it -> [ "${it[1]}.${it[2]}.test.tsv", it.last() ] }

    emit:
    BURDENTEST.out
}

workflow {
    // Define input channels
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.type, row.size, 
            row.category,
            file(row.file)
        ] }

    test_gene_burden(cohorts_ch, model_ch)
}
