#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PENNCNV }   from '../modules/penncnv/penncnv.nf'
include { QUANTISNP } from '../modules/quantisnp/quantisnp.nf'
include { RGADA }     from '../modules/rgada/rgada.nf'
include { CONVERT }   from '../modules/penncnv/convert.nf'
include { SPLIT }     from '../modules/rocker/split.nf'

workflow call_alternates {
    take: 
    signal
    pfb

    main:
    // Prepare signal channel
    tools = Channel.of(params.tools.split(','))
    signal
        | combine(tools)
        | branch {
            merged     : it[2] == 'merged'
            not_merged : params.adjust ? it[2] == 'adjusted' : it[2] == 'raw'
        }
        | set { signal }

    // PENNCNV
    signal.not_merged
        | filter { it.last() == 'penncnv' }
        | combine(pfb)
        | PENNCNV
        | set { penncnv }

    // QUANTISNP
    signal.merged
        | filter { it.last() == 'quantisnp' }
        | QUANTISNP
        | set { quantisnp }

    // RGADA
    signal.merged
        | filter { it.last() == 'rgada' }
        | RGADA
        | set { rgada }

    // Combine calls
    penncnv
        | concat(quantisnp)
        | concat(rgada)
        | map { cohort, key, tool, type, file, nmarker -> [ cohort, key, tool, type.split(',').toList(), file, nmarker ]}
        | transpose
        | collectFile(storeDir: "${params.output_dir}/") { [ "${it[2]}/${it[0]}.${it[2]}.${it[3]}", it[4]] }
        | map { file -> [ file.name.split('\\.')[0], file.name.split('\\.')[1], file.name.split('\\.')[2], file, file.countLines() ] }
        | filter { it.last().toInteger() > 0 }
        | branch {
            // CNV calls
            penn_cnv    : it[2] == 'cnv' && it[1] == 'penncnv'
            notpenn_cnv : it[2] == 'cnv' && it[1] != 'penncnv'

            // Other
            gn  : it[2] == 'gn'
            log : it[2] == 'log'
            qc  : it[2] == 'qc'
        }
        | set { combined_calls }

    // Convert into PennCNV format
    combined_calls.notpenn_cnv
        | combine(pfb)
        | CONVERT
        | set { converted_calls }
    
    // Split
    combined_calls.penn_cnv
        | concat(converted_calls)
        | SPLIT
        | map { cohort, tool, type, cnv, log, nmarker -> [ cohort, tool, type.split(',').toList(), cnv, log, nmarker.split(',').toList() ]}
        | transpose
        | filter { it.last().toInteger() > 1 }
        | set { cnv_loh }

    // Create a combined genotype file for (all) cohorts
    combined_calls.gn
        | collectFile(storeDir: "${params.output_dir}/") { [ "${it[1]}/all.${it[1]}.${it[2]}", it[3]] }
        | map { file -> [ file.name.split('\\.')[0], file.name.split('\\.')[1], file.name.split('\\.')[2], file, file.countLines() ] }
        | concat(combined_calls.gn)
        | set { genotypes }

    emit:
    calls     = cnv_loh
    genotypes = genotypes
    logs      = combined_calls.log
    reports   = combined_calls.qc
}

workflow {
    signal_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }

    pfb_ch    = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }

    call_alternates(signal_ch, pfb_ch)
}
