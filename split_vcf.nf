#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { EXTRACT }     from '../modules/bcftools/extract.nf'
include { CONVERT }     from '../modules/bcftools/convert.nf'
include { SORT_UNIQ }   from '../modules/utils/sort_uniq.nf'

workflow split_vcf {
    take:
    cohort_info

    main:
    cohort_info
        | EXTRACT
        | map { it.last() }
        | collectFile(
            storeDir : "${params.output_dir}/variants",
            name     : 'variants.txt'
        )
        | SORT_UNIQ
        | splitText(
            by   : params.chunk.toInteger(),
            limit: params.limit.toInteger(),
            file : 'chunk',
            elem : 1
        )
        | map { [it.fileName, it] }
        | CONVERT

    emit:
    CONVERT.out
}

// Workflow
workflow {
    // Define input from file
    cohort_info_ch = Channel.fromPath(params.cohort_info)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.file), file(row.index) ]}

    split_vcf(cohort_info_ch)
}
