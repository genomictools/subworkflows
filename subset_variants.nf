#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SELECT }      from '../modules/select.nf'
include { SUBSET }      from '../modules/subset.nf'

chroms_ch   = Channel.of (1..22) | map { "chr$it" }

workflow subset_variants {
    take:
    dbsnp
    cohorts

    main:
    dbsnp
        | combine(chroms_ch)
        | SELECT
        | map { it.last() }
        | splitText(by: params.chunk, file: true) \
        | map { it -> 
            // Split the file name and extract the first par
            def chrom = it.baseName.split("\\.")[1]
            def chunk = it.baseName.split("\\.")[2]
            return [chrom, chunk, it]
        }
        | combine(cohorts)
        | SUBSET
        | filter { it.last().toInteger() > 0 }
        | set { variants }

    emit:
    variants = variants
}

// worflow
workflow {
    dbsnp_ch = Channel.fromFilePairs(params.dbsnp, flat: true) 
        | map { ['dbsnp', it[1], it[2]] }

    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ 
            row.cohort, row.type, row.size,
            file(row.vars_file), file(row.vars_index),
            file(row.population)
        ] }

    subset_variants(dbsnp_ch, cohorts_ch)
}
