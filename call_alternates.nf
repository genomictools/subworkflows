#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PENNCNV }   from '../modules/penncnv/penncnv.nf'
include { QUANTISNP } from '../modules/quantisnp/quantisnp.nf'
include { RGADA }     from '../modules/rgada/rgada.nf'
include { CONVERT }   from '../modules/penncnv/convert.nf'
include { COMBINE }   from '../modules/rocker/combine.nf'
include { PLINK }     from '../modules/plink/plink.nf'

workflow call_alternates {
    take: 
    signal
    genotype
    pfb

    main:
    types = Channel.of(params.type.split(','))
    tools = Channel.of(params.tools.split(','))

    // PLINK
    genotype
        | combine(tools) | filter { it.last() == 'plink' }
        | combine(types) | filter { it.last() == 'roh' }
        | PLINK
        | set { plink }

    // PENNCNV
    signal
        | ( params.adjust ? filter { it[2] == 'adjusted' } : filter { it[2] == 'raw' } )
        | combine(tools) | filter { it.last() == 'penncnv' }
        | combine(types) | filter { it.last() != 'roh' }
        | combine(pfb)
        | PENNCNV
        | set { penncnv }

    // QUANTISNP
    signal
        | filter { it[2] == 'merged' }
        | combine(tools) | filter { it.last() == 'quantisnp' }
        | combine(Channel.of("cnv,loh"))
        | combine(Channel.value(file(params.quant_params)))
        | combine(Channel.value(file(params.quant_ratios)))
        | QUANTISNP
        | set { quantisnp }

    // RGADA
    signal
        | filter { it[2] == 'merged' }
        | combine(tools) | filter { it.last() == 'rgada' }
        | combine(Channel.of("cnv,loh"))
        | RGADA
        | set { rgada }

    // Combine calls
    quantisnp
        | concat(rgada)
        | map {[
            it[0], it[1], it[2],
            it[3].split(',').toList(),
            it[4], it[5],
            it[6].split(',').toList()
        ]}
        | transpose
        | combine(pfb)
        | CONVERT
        | concat(penncnv)
        | groupTuple(by: [0, 2, 3])
        | COMBINE
        | filter { it.last().toInteger() > 1 }
        | set { calls }

    emit:
    calls
}

workflow {
    signal_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }

    genotype_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.level, file(row.file), file(row.log), row.nmarkers ] }

    pfb_ch      = Channel.fromPath(params.pfb) | map { [ it.simpleName, it ] }

    call_alternates(signal_ch, genotype_ch, pfb_ch)
}
