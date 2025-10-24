#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SUBSETMULTIVCF as SUBSET }  from '../modules/bcftools/subset.nf'
include { ANNOTATE }    from '../modules/bcftools/annotate.nf'
include { COMBINE }     from '../modules/bcftools/combine.nf'

workflow get_cohort {
    take:
    cohorts
    bed
    annotations

    main:
    cohorts
        | combine(bed, by: 0)
        | SUBSET        
        | filter { it.last().toInteger() > 0 }
        | ( params.annotate ? combine(annotations, by: [0, 1]) : map { it } )
        | ( params.annotate ? ANNOTATE : map { it } )
        | groupTuple(by: 0)
        | COMBINE
        | filter { it.last().toInteger() > 0 }
        | set { variants }

    emit:
    variants
}

workflow  {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort,
            file(row.file), file(row.index),
            file(row.pedigree)
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
        
    annotations_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ cohort: row.cohort, chrom: row.chrom, start: row.start, end: row.end, file: file(row.annot_file), index: file(row.annot_index) ] }
        | map { it -> 
            chrom = it.chrom ?: (1..2).collect { "chr$it" } + ['chrX', 'chrY']
            key   = (it.start && it.end) ? "${chrom}:${it.start}-${it.end}" : chrom
            [ it.cohort, key, it.file, it.index ]
        }
        | transpose
        | unique

    get_cohort( cohorts_ch, coords_ch, annotations_ch )
}
