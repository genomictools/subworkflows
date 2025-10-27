#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { DOWNLOAD }    from '../modules/gdc/download.nf'
include { MERGE }       from '../modules/bcftools/merge.nf'
include { CONVERT }     from '../modules/pysam/convert.nf'

workflow download_stjude {
    take:
    manifest
    samplesheet
    
    main:
    manifest_ch = Channel.fromPath(manifest)
        | splitCsv(header: false, sep: '\t', skip: 1)
        | map { row -> [
            row[params.manifest_id.toInteger()  - 1],
            row[params.manifest_url.toInteger() - 1]
        ]}

    samplesheet_ch = Channel.fromPath(samplesheet)
        | splitCsv(header: false, sep: '\t', skip: 1)
        | map { row -> [
            row[params.sheet_name.toInteger()  - 1],
            params.cohort
        ]}
        | unique()

    manifest_ch | DOWNLOAD | set { download_ch }

    if ( params.type == 'gvcf') {
    download_ch
        | groupTuple(by: 0, sort: true)
        | CONVERT
        | combine(samplesheet_ch, by: 0)
        | groupTuple(by: 3)
        | MERGE
    }

    emit:
    download_ch
}

workflow {
    download_stjude(params.manifest, params.samplesheet) 
}
