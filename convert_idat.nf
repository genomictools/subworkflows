#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { IDAT2GTC }   from '../modules/bcftools/idat2gtc.nf'
include { GTC2VCF }    from '../modules/bcftools/gtc2vcf.nf'
include { VCF2TXT }    from '../modules/bcftools/vcf2txt.nf'

workflow convert_idat {
    take:
    cohorts

    main:
    cohorts
        | IDAT2GTC
        | GTC2VCF
        | VCF2TXT
        | set { txt }

    emit:
    txt
}

workflow  {
    // Define input channels
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }

    convert_idat(cohorts_ch)
}
