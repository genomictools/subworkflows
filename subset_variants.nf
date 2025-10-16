#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIM }        from '../modules/bcftools/trim.nf'
include { CONVERT }     from '../modules/plink/convert.nf'
include { FILTER }      from '../modules/plink/filter.nf'
include { REMOVE }      from '../modules/plink/remove.nf'
include { COMBINE }     from '../modules/plink/combine.nf'

test_ch = Channel.of(params.tests.split(','))

workflow subset_variants {
    take:
    genotypes
    pedigree

    main:
    genotypes
        | ( params.trim ? TRIM : map { it } )
        | filter { it.last().toInteger() > 0 }
        | combine(pedigree, by: 0)
        | CONVERT
        | filter { it.last().toInteger() > 0 }
        | ( params.filter ? FILTER : map { it } )
        | filter { it.last().toInteger() > 0 }
        | groupTuple(by: [0, 2])
        | COMBINE
        | ( params.remove  ? REMOVE  : map { it } )
        | filter { it.last().toInteger() > 0 }
        | set { combined }

    emit:
    combined
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

    subset_variants( genotypes_ch, pedigree_ch )
}