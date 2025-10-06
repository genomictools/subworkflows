#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { DOWNLOAD }    from '../modules/download.nf'
include { MERGE }       from '../modules/merge.nf'
include { CONVERT }     from '../modules/convert.nf'

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

    manifest_ch
        | DOWNLOAD
        | groupTuple(by: 0, sort: 'hash')
        | CONVERT
        | combine(samplesheet_ch, by: 0)
        | groupTuple(by: 2)
        | MERGE

    emit:
    MERGE.out
}

workflow {
    download_stjude(params.manifest, params.samplesheet) 
}
