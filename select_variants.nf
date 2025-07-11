#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SUBSET }      from '../modules/subset.nf'
include { FILTER }      from '../modules/filter.nf'

workflow select_variants {
    take:
    cohorts
    
    main:
    cohorts
        | SUBSET
        | FILTER

    emit:
    genotypes = FILTER.out
}

workflow  {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, 
            (row.start && row.end) ? "${row.chrom}:${row.start}-${row.end}" : row.chrom,
            file(row.file), file(row.index), file(row.samples)
        ] }

    select_variants(cohorts_ch)
}
