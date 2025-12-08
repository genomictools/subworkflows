#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/rocker/extract.nf'
include { MERGE }     from '../modules/rocker/merge.nf'
include { ADJUST }    from '../modules/penncnv/adjust.nf'
include { MAPSIGNAL } from '../modules/penncnv/mapsignal.nf'
include { CONVERT }   from '../modules/samtools/convert.nf'

workflow prepare_signal {
    take: 
    signal
    pfb
    gcm
    // fasta

    main:
    // Initialize channels
    raw = adjust = merged = Channel.empty()

    // Map & Branch signal by type
    signal
        | branch { 
            text   : it[2] == 'raw' || it[2] == 'gtc' || it[2] == 'report'
            binary : it[2] == 'bam' || it[2] == 'cram'
        }
        | set { signal }

    // Text
    signal.text
        | map { cohort, key, level, file ->
            def log = Path.of("${cohort}.${key}.raw.log")
            def nmarkers = file.readLines().size()
            [ cohort, key, level, file, log, nmarkers ]
        }
        | branch { 
            raw  : it[2] == 'raw'
            gtc  : it[2] == 'gtc'
        }
        | set { text }

    // Process text files
    text.gtc
        | EXTRACT
        | concat(text.raw)
        | ( params.adjust ? combine(gcm) : identity() )
        | ( params.adjust ? ADJUST       : identity() )
        | ( params.tools.contains('quantisnp') || params.tools.contains('rgada') ? combine(pfb) : identity() )
        | ( params.tools.contains('quantisnp') || params.tools.contains('rgada') ? MERGE : identity() )
        | filter { it.last().toInteger() > 1 }
        | set { text }

    Channel.empty()
        | ( params.adjust ? concat(ADJUST.out) : identity() )
        | ( params.tools.contains('quantisnp') || params.tools.contains('rgada') ? concat(MERGE.out) : identity() )
        | set { text }

    // Binary
    intervals = Channel.fromPath(params.intervals)
        | splitCsv(header: false, sep: '\t')
        | map { chrom, start, end -> [ chrom  + ":" + (start.toInteger() +1) + "-" + (end.toInteger() +1) ] }

    signal.binary
        | map { cohort, key, level, file ->
            def log = Path.of("${cohort}.${key}.raw.log")
            [ cohort, key, level, file.sort(), log, 0 ]
        }
        | ( params.intervals != null ? combine(intervals) : identity() )
        | branch { 
            bam  : it[2] == 'bam'
            cram : it[2] == 'cram'
        }
        | set { binary }

    binary.cram
        | CONVERT
        | concat(binary.bam)
        | combine(pfb)
        | MAPSIGNAL
        | collectFile(storeDir: "${params.output_dir}/mappedsignal", keepHeader: true) { [ "${it[0]}.${it[1]}.${it[2]}.txt", it[3]] }
        | map { file ->
            def cohort = file.name.split('\\.')[0]
            def key = file.name.split('\\.')[1]
            def level = file.name.split('\\.')[2]
            def log = Path.of("${file.name}.log")
            def nmarkers = file.countLines()
            [ cohort, key, level, file, log, nmarkers ] 
        }
        | view
        | set { binary }

    // Combine
    text
        // | concat(binary)
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