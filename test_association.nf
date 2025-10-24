#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ASSOCTEST }       from '../modules/plink/assoctest.nf'
include { PLOTMANHATTAN }   from '../modules/rocker/plotmanhattan.nf'

workflow test_association {
    take:
    genotypes
    phenotypes
    covariates
    tests

    main:
    genotypes
        | combine(tests)
        | combine(phenotypes, by: 0)
        | combine(covariates, by: [0,1])
        | ASSOCTEST
        | transpose
        | map { it -> 
            def phenotype = it[3].name.split('\\.')[2]
            [ it[0], it[1], it[2], phenotype, it[3], it[4] ]
        }
        | ( params.plot ? PLOTMANHATTAN : map { it })

    emit:
    tests = ASSOCTEST.out
}

workflow  {
    genotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.category,
            file(row.bim),file(row.bed), file(row.fam),file(row.log),
            row.n_samples, row.n_variants
        ] }

    phenotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort,file(row.phenotypes) ] }
        | unique

    covariates_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.category, file(row.covariates) ] }
        | unique

    test_ch = Channel.of(params.tests.split(','))

    test_association( genotypes_ch, phenotypes_ch, covariates_ch, test_ch )
}