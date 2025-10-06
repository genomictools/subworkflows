#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT }     from '../modules/plink/convert.nf'
include { EXTRACT }     from '../modules/plink/extract.nf'
include { AGGREGATE }   from '../modules/rocker/aggregate.nf'

category_ch = Channel.of(params.categories.split(','))
variable_ch = Channel.of(params.variables.split(','))

workflow summarize_cases {
    take:
    genotypes
    pedigree
    annotations

    main:
    genotypes
        | combine(pedigree, by: 0)
        | CONVERT
        | combine(variable_ch)
        | EXTRACT
        | filter { it[3] == 'rlist' }
        | join(annotations, by: [0,1,2])
        | AGGREGATE
        | concat(EXTRACT.out)
        | concat(annotations)
        | set { summary }

    emit:
    summary = summary
}

workflow  {
    genotypes_ch = Channel.fromPath(params.genotypes)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.category, row.file, row.index, row.n_vars ] }
    
    pedigree_ch = Channel.fromPath(params.pedigree)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.phenotype) ] }
        | unique

    annotations_ch = Channel.fromPath(params.annotations)
        | map { row -> [ row.cohort, row.key, row.category, row.variable, row.file ] }

    summarize_cases( genotypes_ch, pedigree_ch, annotations_ch )
}
