#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from '../modules/coordinates.nf'

workflow get_coordinates {
    take:
    coords
    genome
    style
    
    main:
    COORDINATES(coords, genome, style)
        | transpose
        | filter { it.last().size() > 0 }
        | set { bed }

        bed
            | collectFile { it -> [ "${it.first()}.bed", it.last() ] } 
            | map { [ it.simpleName, it ] }
            | splitText(
                by: (params.coding ? params.chunk.toInteger() : 1),
                file: 'chunk'
            )
            | map { [ it.first(), it.last().fileName, it.last() ] }
            | set { chunks }

    emit:
    bed    = bed
    chunks = chunks
}

workflow  {
    coords_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ cohort: row.cohort, chrom: row.chrom, start: row.start, end: row.end, genelist: file(row.genelist) ] }
        | map { it -> 
            chrom = it.chrom ?: (1..22).collect { "chr$it" } + ['chrX', 'chrY']
            key   = (it.start && it.end) ? "${chrom}:${it.start}-${it.end}" : chrom
            [ it.cohort, key, chrom, it.start ?: null, it.end ?: null, it.genelist ? file(it.genelist) : null ]
        }
        | transpose
        | unique
        | groupTuple(by: [1,2,3,4])
    
    get_coordinates( coords_ch, params.genome, params.style )
}
