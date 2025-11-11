#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { DOWNLOAD }    from '../modules/dx/download.nf'

workflow download_stjude {
    take:
    manifest

    main:
    manifest
        | DOWNLOAD
        | set { download_ch }

    emit:
    download_ch
}

workflow {
    manifest_ch = Channel.fromPath(params.manifest)
        | splitCsv(header: true)
        | map { row -> [ row.project, row.directory, row.filename ] }

    download_stjude( manifest_ch )
}
