#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT }   from '../modules/rocker/extract.nf'
include { MERGE }     from '../modules/rocker/merge.nf'
include { ADJUST }    from '../modules/penncnv/adjust.nf'
include { MAPSIGNAL } from '../modules/penncnv/mapsignal.nf'
include { CONVERT }   from '../modules/samtools/convert.nf'

include { IDAT2GTC }   from '../modules/bcftools/idat2gtc.nf'
include { GTC2VCF }    from '../modules/bcftools/gtc2vcf.nf'
include { VCF2TXT }    from '../modules/bcftools/vcf2txt.nf'

workflow prepare_signal {
    take: 
    signal
    pfb
    gcm

    main:
    // Initialize channels
    raw = adjust = merged = Channel.empty()

    // Map & Branch signal by type
    def file_types = branchCriteria { it ->
        raw     : it[2] == 'raw'
        report  : it[2] == 'report'

        idat : it[2] == 'idat'
        gtc  : it[2] == 'gtc'
        vcf  : it[2] == 'vcf'

        bam  : it[2] == 'bam'
        cram : it[2] == 'cram'
    }

    signal
        | map { cohort, key, level, file ->
            def files = file instanceof List ? file.sort() : file
            def log = Path.of("${cohort}.${key}.raw.log")
            def nmarkers = ['raw', 'report'].contains(level) ? file.readLines().size() : 0
            [ cohort, key, level, files, log, nmarkers ]
        }
        | branch( file_types )
        | set { signal }

    // Report
    signal.report
        | EXTRACT
        | set { report }

    // Illumina
    signal.idat
        | IDAT2GTC
        | concat(signal.gtc)
        | GTC2VCF
        | concat(signal.vcf)
        | VCF2TXT
        | concat(report)
        | concat(signal.raw)
        | ( params.adjust ? combine(gcm) : identity() )
        | ( params.adjust ? ADJUST       : identity() )
        | set { text }

    text
        | ( params.tools.contains('quantisnp') || params.tools.contains('rgada') ? combine(pfb) : identity() )
        | ( params.tools.contains('quantisnp') || params.tools.contains('rgada') ? MERGE : identity() )
        | set { merged }

    // Binary
    intervals = Channel.fromPath(params.intervals)
        | splitCsv(header: false, sep: '\t')
        | map { chrom, start, end -> [ chrom  + ":" + (start.toInteger() +1) + "-" + (end.toInteger() +1) ] }

    signal.bam
        | concat(signal.cram)
        | ( params.intervals != null ? combine(intervals) : identity() )
        | branch( file_types )
        | set { binary }

    binary.cram
        | CONVERT
        | concat(binary.bam)
        | combine(pfb)
        | MAPSIGNAL
        | collectFile(storeDir: "${params.output_dir}/mappedsignal", keepHeader: true) { [ "${it[0]}.${it[1]}.${it[2]}.txt", it[3]] }
        | map { file ->
            def file_name = file.name.split('\\.')
            def log = Path.of("${file.name}.log")
            def nmarkers = file.countLines()
            [ file_name[0], file_name[1], file_name[2], file, log, nmarkers ] 
        }
        | set { binary }

    // Combine
    text
        | concat(merged)
        | concat(binary)
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