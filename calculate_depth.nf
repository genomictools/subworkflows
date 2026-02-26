#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from '../modules/bioconductor/coordinates.nf'
include { DEPTH }       from '../modules/bamsignals/depth.nf'
include { STATS }       from '../modules/samtools/stats.nf'

workflow calculate_depth {
    take:
    cohorts
    genes

    main:
    if ( params.coords == null ) {
    genes
        | concat(Channel.of("reference"))
        | combine(Channel.of("exon", "gene"))
        | COORDINATES
        | set { coordinates }
    } else {
    Channel.fromPath(params.coords)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.gene, row.feature, file(row.coords) ] }
        | set { coordinates }
    }

    cohorts
        | combine(coordinates)
        | DEPTH
        | branch { 
            coverage : it[4] == 'gene'
            counts   : it[4] == 'exon'
        }
        | set { depth }

    cohorts
        | STATS
        | transpose
        | map { cohort, key, type, stats -> [ cohort, key, type, stats.name.split('\\.')[3], stats ] }
        | filter { it[3] == 'flagstats' }
        | set { stats }

    emit:
        coverage = depth.coverage
        counts   = depth.counts
        stats    = stats
}

workflow {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, file(row.bam), file(row.bai) ] }

    genes_ch = Channel.fromPath(params.genes)
        | splitText
        | map { it.trim() }

    calculate_depth(cohorts_ch, genes_ch )
}