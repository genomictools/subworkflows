#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'
include { MERGE }     from '../modules/merge.nf'
include { GENOTYPE }  from '../modules/genotype.nf'

workflow prepare_signal {
    take: 
    gtc
    pfb
    gcm

    main:
    // Extract genotype
    genotype = Channel.empty()
    if ( params.genotype ) {
    gtc
        | combine(pfb)
        | GENOTYPE
        | set { genotype }
    }

    // Extract (and adjust) signal
    gtc
        | ( params.extract ? EXTRACT : map {
            cohort, key, file ->
            def nmarkers = new File(file.toString()).readLines().size()
            def log      = Path.of(file.simpleName + ".log")
            [ cohort, key, 'raw', file, log, nmarkers ]
        } )
        | filter { it.last().toInteger() > 1 }
        | set { raw }

    adjust = Channel.empty()
    if ( params.adjust ) {
    raw
        | combine(gcm)
        | ADJUST
        | filter { it.last().toInteger() > 1 }
        | set { adjust }
    }

    merge = Channel.empty()
    if ( params.merge ) {
    raw
        | combine(pfb)
        | MERGE
        | filter { it.last().toInteger() > 1 }
        | set { merge }
    }

    // Combine
    raw
        | concat(adjust)
        | concat(merge)
        | set { signal }

    emit:
    signal
    genotype
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, file(row.file) ] }
    pfb    = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm    = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }

    prepare_signal(gtc_ch, pfb, gcm)
}