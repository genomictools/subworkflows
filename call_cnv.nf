#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXONS }   from '../modules/bioconductor/exons.nf'
include { COUNT }   from '../modules/exomedepth/count.nf'
include { CALL }    from '../modules/exomedepth/call.nf'
include { PLOT }    from '../modules/exomedepth/plot.nf'
include { REPORT }  from '../modules/rocker/report.nf'
include { VISUALIZE } from '../modules/bamsignals/visualize.nf'

workflow call_cnv {
    take:
    cohorts
    stats
    genes

    main:
    // Extract exons for the given gene (exonic) ranges
    genes
        | EXONS
        | set { exons }

    // Count reads in exonic regions for all samples
    exons
        | combine( cohorts )
        | combine( Channel.fromPath(params.fasta) )
        | COUNT
        | filter { it.last().toInteger() > params.coverage }
        | branch {
            cases    : it[2] == 'case'
            controls : it[2] == 'control'
        }
        | set { counts }

    // Aggregate reference ranges, in each case sample
    stats
        | map { [ it.cohort, it.key, it['ref.samples']] }
        | combine( counts.cases, by: [0,1] )    
        | map { cohort, key, ref_key, type, gene, count, coverage ->
            [ cohort, ref_key, type, gene, count, coverage, key ]
        }
        | combine( counts.controls, by: [0,1,3] )
        | set { combined_counts }
    
    // Call CNVs
    combined_counts
        | groupTuple(by: [0,2,4,6])
        | map { cohort, ref_key, gene, type, counts, coverage, key, ref_type, ref_counts, ref_coverage ->
            [ cohort, gene, key, counts, ref_counts ]
        }
        | CALL
        | filter { it.last().toInteger() > 1 }
        | set { calls }

    // Save output to file
    calls
        | collectFile(
            storeDir: "${params.output_dir}/calls",
            keepHeader: true,
            skip: 1
        ) { [ "${it[0]}.cnv.tsv", it[4]] }
        | map { file -> [ file.name.split('\\.')[0], 'calls', file, file.readLines().size()] }
        | filter { it.last() > 1 }
        | REPORT
        | set { reports }

    // Collect output in a file & filter
    calls
        | map { it[0,1,2,4]}
        | splitCsv(header: true, sep: '\t')
        | map { cohort, gene, key, cnv -> [ cohort: cohort, gene: gene, key: key ] + cnv }
        | filter { it['BF'].toDouble() >= params.BF }
        | filter { it['reads.ratio'] != 'NA' && it['reads.ratio'] != 'Inf'}
        | filter { it['reads.ratio'].toDouble() >= params.dup_ratio || it['reads.ratio'].toDouble() <= params.del_ratio}
        | map { [it.cohort, it.gene, it.key] }
        | set { filtered_calls }

    // Plot genes with CNVs
    if ( params.plot_counts ) {
        calls
            | combine( filtered_calls, by: [0,1,2] )
            | PLOT
    }

    // Visualize bams
    if ( params.plot_coverage ) {
        combined_counts
        | map { cohort, ref_key, gene, type, counts, coverage, key, ref_type, ref_counts, ref_coverage -> [ cohort, gene, key, ref_key ] }
        | combine( filtered_calls, by: [0,1,2] )
        | map { cohort, gene, key, ref_key -> [ gene, cohort, key, ref_key ] }
        | combine(exons, by: 0)
        | map { gene, cohort, key, ref_key, exons -> [ cohort, [key, ref_key], key, gene, exons ] }
        | transpose
        | combine(cohorts, by: [0,1])
        | unique
        | groupTuple(by: [0,3,2,4])
        | map { cohort, all_keys, key, gene, exons, type, bam, bai ->[ cohort, gene, key, exons, bam, bai ]}
        | VISUALIZE
    }

    emit:
    calls   = calls
    reports = reports
}

workflow {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, file(row.bam), file(row.bai) ] }
    
    ref_ch = Channel.fromPath(params.refs)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.ref ] }

    genes_ch = Channel.fromPath(params.genes)
        | splitText
        | map { it.trim() }

    call_cnv(cohorts_ch, stats_ch, genes_ch )
}