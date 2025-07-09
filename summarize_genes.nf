#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CONVERT }     from '../modules/convert.nf'
include { EXTRACT }     from '../modules/extract.nf'
include { AGGREGATE }   from '../modules/aggregate.nf'

variable_ch = Channel.of( 'rlist', 'snplist', 'frqx', 'frq.strat' )

workflow summarize_genes {
    take:
    genotypes
    phenotypes
    annotations

    main:
    genotypes
        | combine(phenotypes, by: 0)
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
    
    phenotypes_ch = Channel.fromPath(params.phenotypes)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.phenotype) ] }
        | unique

    annotations_ch = Channel.fromPath(params.annotations)
        | map { row -> [ row.cohort, row.key, row.category, row.variable, row.file ] }

    summarize_genes( genotypes_ch, phenotypes_ch, annotations_ch )
}
