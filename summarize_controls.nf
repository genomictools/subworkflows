#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SUBSET }      from '../modules/bcftools/subset.nf'
include { FILTER }      from '../modules/bcftools/filter.nf'
include { TABULATE }    from '../modules/bcftools/tabulate.nf'
include { SUMMARIZE }   from '../modules/rocker/summarize.nf'

workflow summarize_controls {
    take:
    cohorts
    categories

    main:
    cohorts
        | ( params.subset ? SUBSET : map { it + [ 1 ] } )
        | filter { it.last().toInteger() > 0 }
        | combine( categories )
        | FILTER
        | filter { it.last().toInteger() > 0 }
        | TABULATE
        | filter { it.last().toInteger() > 0 }
        | groupTuple( by: [0,2] )
        | SUMMARIZE
        | set { summary }

    emit:
    summary = summary
}

workflow  {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.key,
            file(row.file), file(row.index)
        ] }

    // 'Pathogenic,Damaging,Splicing,High,PTV,Stop,Rare'
    category_ch = Channel.of(params.categories.split(','))
    
    summarize_controls( 
        cohorts_ch,
        category_ch
    )
}
