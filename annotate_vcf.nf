#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { ANNOTATE }    from '../modules/annotate.nf'
include { FORMAT }      from '../modules/format.nf'

workflow annotate_vcf {
    take:
    cohort_info
    annotations

    main:
    annotations
        | combine(cohort_info)
        | ANNOTATE
        | FORMAT

    emit:
    ANNOTATE.out
}

// Workflow
workflow {
    // Define input from file
    cohort_info_ch = Channel.fromPath(params.cohort_info)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.file), file(row.index) ]}

    annotations_ch = Channel.fromPath(params.annotations)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.assembly, row.tool, row.version,
            file(row.file), file(row.index)
        ]}

    annotate_vcf(cohort_info_ch, annotations_ch)
}
