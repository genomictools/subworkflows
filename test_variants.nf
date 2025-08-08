#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT }     from '../modules/convert.nf'
include { PRUNE }       from '../modules/prune.nf'
include { COMBINE }     from '../modules/combine.nf'
include { FILTER }      from '../modules/filter.nf'
include { TEST }        from '../modules/test.nf'
include { PLOT }        from '../modules/plot.nf'

test_ch = Channel.of(params.tests.split(','))

workflow test_variants {
    take:
    genotypes
    pedigree
    phenotypes

    main:
    genotypes
        | combine(pedigree, by: 0)
        | CONVERT
        | ( params.prune ? PRUNE : map { it } )
        | groupTuple(by: [0, 2])
        | COMBINE
        | ( params.filter ? FILTER : map { it } )
        | combine(test_ch)
        | combine(phenotypes, by: 0)
        | TEST
        | transpose
        | map { it -> 
            def phenotype = it[3].name.split('\\.')[2]
            [ it[0], it[1], it[2], phenotype, it[3], it[4], it[5] ]
        }
        | ( params.plot ? PLOT : map { it })

    emit:
    tests = TEST.out
}

workflow  {
    genotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort,row.key,row.category,file(row.file),file(row.index),
            row.n_samples,row.n_variants
        ] }

    pedigree_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort,file(row.pedigree) ] }
        | unique

    phenotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort,file(row.phenotypes) ] }
        | unique

    test_variants( genotypes_ch, pedigree_ch, phenotypes_ch )
}
