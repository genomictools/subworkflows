#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MERGE }       from '../modules/plink/merge.nf'
include { FILTER }      from '../modules/plink/filter.nf'
include { SAMPLE }      from '../modules/plink/sample.nf'
include { SCALE }       from '../modules/plink/scale.nf'
include { ASSIGN }      from '../modules/rocker/assign.nf'
include { PLOTPCA }     from '../modules/rocker/plotpca.nf'

modes_ch    = Channel.of(params.modes.split(','))

workflow infer_ancestry {
    take:
    cases
    references
    population

    main:
    marged = MERGE(cases, references)
    marged
        | FILTER
        | SAMPLE
        | combine(population)
        | combine(modes_ch)
        | SCALE
        | ASSIGN
        | PLOTPCA

    emit:
    plots = PLOTPCA.out
}

// worflow
workflow {
    cases_ch = Channel.fromPath(params.cases)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.type, file(row.bim), file(row.bed), file(row.fam), file(row.nosex), file(row.log)] }
    
    reference_ch = Channel.fromPath(params.cases)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.type, file(row.bim), file(row.bed), file(row.fam), file(row.nosex), file(row.log)] }

    population_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.population) ] }

    infer_ancestry(cases_ch, reference_ch, population_ch)
}
