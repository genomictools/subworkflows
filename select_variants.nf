#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SPLIT }       from '../modules/split.nf'
include { FILL }        from '../modules/fill.nf'
include { FILTER }      from '../modules/filter.nf'

category_ch = Channel.of(params.categories.split(','))

workflow select_variants {
    take:
    variants
    chunks

    main:
    variants
        | combine(chunks, by: 0)
        | SPLIT
        | ( params.fill ? FILL : map {it} )
        | filter { it.last().toInteger() > 0 }
        | combine(category_ch)
        | FILTER
        | filter { it.last().toInteger() > 0 }
        | multiMap {
            genotypes   : [ it[0], it[1], it[2], it[3], it[4], it[7], it[8] ]
            annotations : [ it[0], it[1], it[2], 'annotations', it[5] ]
            qc          : [ it[0], it[1], it[2], 'qc', it[6] ]
        }
        | set { selected }

    emit:
    genotypes   = selected.genotypes
    annotations = selected.annotations
    qc          = selected.qc
}

workflow  {
    variants_ch = Channel.fromPath(params.variants)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.file, row.index, row.n_samples, row.n_variants ] }

    chunks_ch = Channel.fromPath(params.chunks)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.file ] }

    select_variants( variants_ch, chunks_ch )
}
