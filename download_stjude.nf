#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { DOWNLOAD }    from '../modules/gdc/download.nf'
include { COMBINE }     from '../modules/gatk/combine.nf'
include { GENOTYPE }    from '../modules/gatk/genotype.nf'
include { GATHER }      from '../modules/gatk/gather.nf'

workflow download_stjude {
    take:
    manifest
    samplesheet
    fasta_ch

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
        | combine(samplesheet_ch, by: 0)
        | map { id, files, cohort -> [ cohort, id ] + files.flatten() }
        | set { gvcf_ch }

    if ( params.combine ) {
    gvcf_ch
        | groupTuple(by: 0, sort: true)
        | combine(fasta_ch)
        | COMBINE
        | set { combined_gvcf_ch }
    }

    if ( params.joint_call ) {
    combined_gvcf_ch
        | combine(fasta_ch, by: [0, 1])
        | GENOTYPE
        | groupTuple(by: [0,2], sort: true)
        | GATHER
        | set { joint_vcf_ch }
    }
    }

    emit:
    download_ch
}

workflow {
    download_stjude(params.manifest, params.samplesheet) 
}
