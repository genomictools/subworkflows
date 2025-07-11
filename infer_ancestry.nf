#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MERGE }       from '../modules/merge.nf'
include { FILTER }      from '../modules/filter.nf'
include { SCALE }       from '../modules/scale.nf'
include { ASSIGN }      from '../modules/assign.nf'
include { PLOT }        from '../modules/plot.nf'

modes_ch    = Channel.of(params.modes.split(','))

workflow infer_ancestry {
    take:
    cases
    references

    main:
    cases
        | combine(references)
        | MERGE
        | FILTER
        | combine(modes_ch)
        | SCALE
        | ASSIGN
        | PLOT
    
    emit:
    plots = PLOT.out
}

// worflow
workflow {
    population_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.population) ] }

    cases_ch = Channel.fromPath(params.cases)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.type, file(row.bim), file(row.bed), file(row.fam), file(row.nosex), file(row.log)] }
    
    reference_ch = Channel.fromPath(params.cases)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.type, file(row.bim), file(row.bed), file(row.fam), file(row.nosex), file(row.log)] }

    infer_ancestry(cases_ch, reference_ch)
}