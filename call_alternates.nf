#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'
include { DETECT }    from '../modules/detect.nf'

type_ch = Channel.of(params.type.split(','))
format_ch = Channel.of( 'bed', 'tab')

workflow call_alternates {
    take: 
    gtc
    pfb
    gcm
    hmm
    hmm0
    genes
    links

    main:
    gtc
        | EXTRACT
        | ( params.adjust ? combine(gcm) : map { it } )
        | ( params.adjust ? ADJUST       : map { it } )
        | combine(pfb)
        | combine(hmm)
        | combine(hmm0)
        | combine(type_ch)
        | DETECT
        | branch {
            cnv : it[2] == 'cnv'
            loh : it[2] == 'loh'
        }
        | set { alternates }

    emit:
    signal = EXTRACT.out
    cnv    = alternates.cnv
    loh    = alternates.loh
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }

    hmm     = Channel.fromPath(params.hmm)
    hmm0    = Channel.fromPath(params.hmm0)
    genes   = Channel.fromPath(params.refgene)
    links   = Channel.fromPath(params.reflink)

    pfb = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }

    call_alternates(gtc_ch, pfb, gcm, hmm, hmm0, genes, links)
}