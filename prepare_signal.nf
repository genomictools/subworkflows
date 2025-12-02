#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/rocker/extract.nf'
include { MERGE }     from '../modules/rocker/merge.nf'
include { ADJUST }    from '../modules/penncnv/adjust.nf'

workflow prepare_signal {
    take: 
    signal
    pfb
    gcm

    main:
    // Initialize channels
    raw = adjust = merged = Channel.empty()

    // Map & Branch signal by level
    signal
        | map { cohort, key, level, file ->
            def log = Path.of("${cohort}.${key}.raw.log")
            def nmarkers = file.readLines().size()
            [ cohort, key, level, file, log, nmarkers ]
        }
        | branch { 
            raw : it[2] == 'raw'
            gtc : it[2] == 'gtc'
        }
        | set { signal }

    // Extract raw signals from gtc
    signal.gtc
        | EXTRACT
        | concat(signal.raw)
        | set { raw }

    // Adjust with GC
    if ( params.adjust ) {
        raw
            | combine(gcm)
            | ADJUST
            | filter { it.last().toInteger() > 1 }
            | set { adjust }
    }

    // Merge with pfb, when tools include quantisnp, or rgada
    if ( params.tools.contains('quantisnp') || params.tools.contains('rgada') ) {
        adjust
            | combine(pfb)
            | MERGE
            | filter { it.last().toInteger() > 1 }
            | set { merged }
    }

    // Combine
    signal.gtc
        | concat(raw)
        | concat(adjust)
        | concat(merged)
        | set { signal }

    emit:
    signal
}

workflow {
    signal_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file) ] }
    pedigree_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.pedigree) ] }
    pfb    = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm    = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }

    prepare_signal(signal_ch, pfb, gcm)
}