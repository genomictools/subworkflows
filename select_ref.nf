#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXONS }   from '../modules/bioconductor/exons.nf'
include { COUNT }   from '../modules/exomedepth/count.nf'
include { SELECT }  from '../modules/exomedepth/select.nf'
include { REPORT }  from '../modules/rocker/report.nf'
include { STATS }   from '../modules/samtools/stats.nf'

workflow select_ref {
    take:
    cohorts

    main:
    // Get stats for all samples
    cohorts
        | STATS
        | transpose
        | map { cohort, key, type, stats ->
                [ cohort, key, type, stats.name.split('\\.')[3], stats ] 
        }
        | set { stats }

    // Extract exons for the given exonic ranges
    Channel.of("reference")
        | EXONS
        | set { exons }
    
    // Count exons in reference samples
    exons
        | combine( cohorts )
        | combine( Channel.fromPath(params.fasta) )
        | COUNT
        | filter { it.last().toInteger() > params.bins * 5}
        | map { cohort, key, type, range, counts, coverage ->
            def new_key = ( type == 'control' ) ? 'control' : key
            [ cohort, new_key, counts ]
        }
        | groupTuple(by: [0,1])
        | branch {
            controls : it[1] == 'control'
            cases    : it[1] != 'control'
        }
        | set { counts }

    // Create a controls set from reference counts
    counts.cases
        | combine(counts.controls, by: 0)
        | SELECT
        | splitCsv(header: true, sep: '\t')
        | map { cohort, key, ref -> [ cohort: cohort, key: key ] + ref }
        | filter { it.selected == 'TRUE' }
        | set { ref }

    // Save output to file
    SELECT.out 
        | collectFile(
            storeDir: "${params.output_dir}/select",
            keepHeader: true,
            skip: 1
        ) { [ "${it[0]}.select.tsv", it[2]] }
        | map { file -> [ file.name.split('\\.')[0], 'samples', file, file.readLines().size()] }
        | filter { it.last() > 1 }
        | REPORT
        | set { reports }

    emit:
    ref     = ref
    stats   = stats
    reports = reports
}

workflow {
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.key, row.type, file(row.bam), file(row.bai) ] }

    select_ref(cohorts_ch)
}