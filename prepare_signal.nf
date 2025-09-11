#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/extract.nf'
include { ADJUST }    from '../modules/adjust.nf'

workflow prepare_signal {
    take: 
    gtc
    pfb
    gcm

    main:
    // Extract signal
    gtc
        | ( params.extract ? combine(pfb) : map { it } )
        | ( params.extract ? EXTRACT : map { it } )
        | map { cohort, key, level, files ->
            def level_list = level.tokenize(',')
            def files_list = files instanceof List ? files : [files]
            def nmark_list = files_list.collect { new File(it.toString()).readLines().size() }
            def log_list   = files_list.collect { Path.of(it.baseName + ".log") }
            [ cohort, key, level_list, files_list, log_list, nmark_list ]
        }
        | transpose
        | filter { it.last().toInteger() > 1 }
        | set { signal }

    // Genotype
    signal
        | filter { it[2] == 'genotype' }
        | set { genotype }

    // Adjust with GC
    adjust = Channel.empty()
    if ( params.adjust ) {
        signal
            | filter { it[2] == 'raw' }
            | combine(gcm)
            | ADJUST
            | filter { it.last().toInteger() > 1 }
            | set { adjust }
    }

    // Combine
    signal
        | concat(adjust)
        | set { signal }

    emit:
    signal
    genotype
}

workflow {
    gtc_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file) ] }
    pfb    = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }
    gcm    = Channel.fromPath(params.gcm) | map { [ it.simpleName, it ] }

    prepare_signal(gtc_ch, pfb, gcm)
}