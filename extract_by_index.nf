#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { EXTRACT }     from '../modules/extract.nf'
include { RETRIEVE }    from '../modules/retrieve.nf'
include { SPLIT }       from '../modules/split.nf'
include { COMBINE }     from '../modules/combine.nf'
include { PAIR }        from '../modules/pair.nf'

// Extract reads from fastq files based on index sequences
workflow extract_by_index {
    take: 
        undetermined
        unkown

    main:
    // Generate ranges channel
    ranges_ch = Channel.from( (1..params.nseq).step(params.chunk).collect { i -> 
            def end = Math.min(i + params.chunk - 1, params.nseq)
            [ i, end ]
        })

    undetermined
        | combine(ranges_ch)
        | SPLIT
        | set { fastq }

    // Get index sequences from fastq files
    unkown
        | splitCsv(header: true, sep: ',')
        | multiMap { cohort, row -> 
            both   : [ cohort, row.Lane, "${row.index}+${row.index2}", row.'# Reads', row.'% of Unknown Barcodes', row.'% of All Reads' ]
            index1 : [ cohort, row.Lane, row.index,                    row.'# Reads',  row.'% of Unknown Barcodes', row.'% of All Reads' ]
            index2 : [ cohort, row.Lane, row.index2,                   row.'# Reads',  row.'% of Unknown Barcodes', row.'% of All Reads' ]
        }
        | set { unknown_index }

    unknown_index.both
        | ( params.index1 ? concat(unknown_index.index1) : map { it } )
        | ( params.index2 ? concat(unknown_index.index2) : map { it } )
        | unique
        | groupTuple(by: [0, 1, 2] )
        | map { [ it[0], [ lane: it[1], index: it[2], nread: it[3].collect { it.toInteger() }.sum()]] }
        | filter { it.last().nread.toInteger() > params.min_read } 
        | combine(fastq, by: 0)
        | EXTRACT
        | filter { it.last().size() > 0 }
        | RETRIEVE
        | filter { it.last().size() > 0 }
        | groupTuple(by: [0, 1, 2, 3])
        | COMBINE
        | groupTuple(by: [0, 1, 2])
        | ( params.paired ? PAIR : map { it } )
        | transpose
        | map { [ it[0], it[1], "${it[4].name.replaceAll(/\.fastq\.gz$/, '')}", it[3], it[4] ] }
        | set { retrieved }

    emit:
        retrieved = retrieved
}

workflow {
    // Get undetermined from cohorts
    undetermined = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, 'Undetermined', row.undetermined.split('/|\\.')[1], row.read, file(row.undetermined) ] }

    // Get unknown barcodes from cohorts
    unknown = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.unknown) ] }

    // Extract reads from fastq files
    extract_by_index(undetermined, unknown)
}