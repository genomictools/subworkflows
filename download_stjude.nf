#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { DOWNLOAD }    from '../modules/gdc/download.nf'

workflow download_stjude {
    take:
    manifest

    main:
    manifest_ch = Channel.fromPath(manifest)
        | splitCsv(header: false, sep: '\t')
        | map { row -> [
            row[params.manifest_id.toInteger()  - 1],
            row[params.manifest_url.toInteger() - 1]
        ]}

    manifest_ch | DOWNLOAD | set { download_ch }

    emit:
    download_ch
}

workflow {
    download_stjude(params.manifest)
}
