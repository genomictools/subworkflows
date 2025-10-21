#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXCLUDE }     from '../modules/plink/exclude.nf'
include { PRUNE }       from '../modules/plink/prune.nf'
include { SCALE }       from '../modules/plink/scale.nf'

workflow scale_variants {
    take:
    genotypes

    main:
    genotypes
        | ( params.prune  ? PRUNE  : map { it } )
        | filter { it.last().toInteger() > 0 }
        | ( params.exclude ? EXCLUDE : map { it } )
        | filter { it.last().toInteger() > 0 }
        | combine(Channel.of("noclusters"))
        | combine(Channel.fromPath(params.exlude_regions))
        | SCALE
        | map { [it[0], it[1], it[4]]}
        | set { scaled }

    emit:
    scaled
}

workflow  {
    genotypes_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.cohort, row.category,
            file(row.file), file(row.index),
            row.n_samples,row.n_variants
        ] }

    scale_variants( genotypes_ch )
}