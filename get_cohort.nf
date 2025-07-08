#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SUBSET }      from '../modules/subset.nf'
include { COMBINE }     from '../modules/combine.nf'

workflow get_cohort {
    take:
    cohorts
    bed

    main:
    cohorts
        | combine(bed, by: 0)
        | SUBSET        
        | filter { it.last().toInteger() > 0 }
        | groupTuple(by: 0)
        | COMBINE
        | filter { it.last().toInteger() > 0 }
        | set { variants}

    emit:
    variants
}

workflow  {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort,
            file(row.file), file(row.index),
            file(row.samples)
        ] }
        | unique
        
    coords_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ cohort: row.cohort, chrom: row.chrom, start: row.start, end: row.end ] }
        | map { it -> 
            chrom = it.chrom ?: (1..2).collect { "chr$it" } + ['chrX', 'chrY']
            key   = (it.start && it.end) ? "${chrom}:${it.start}-${it.end}" : chrom
            [ it.cohort, key, chrom, it.start ?: null, it.end ?: null]
        }
        | transpose
        | unique
        | groupTuple(by: [1,2,3,4])

    get_cohort( cohorts_ch, coords_ch )
}
