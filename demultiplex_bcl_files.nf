#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { BCLCONVERT }  from '../modules/bclconvert.nf'

// Demultiplexing bcl files
workflow demultiplex_bcl_files {
    take: 
        cohorts
    
    main:
    // Convert bcl files to fastq
    cohorts
        | BCLCONVERT
        | multiMap {
            fastq   : [ it[0], it[1] ]
            reports : [ it[0], it[2] ]
            logs    : [ it[0], it[3] ]
        }
        | set { bcl_out }
    
    // Split fastq files into determined and undetermined
    bcl_out.fastq
        | transpose
        | map { [
            it.first(),
            it.last().simpleName.split('_')[0] == 'Undetermined' ? 'Undetermined' : 'Determined',
            it.last().simpleName,
            it.last().simpleName.split('_')[3],
            it.last() 
        ] }
        | branch {
            undetermined : it[1] == 'Undetermined'
            determined   : it[1] != 'Undetermined'
        }
        | set { fastq }
    
    bcl_out.reports
        | transpose 
        | branch { 
            known:   it.last().simpleName == 'Demultiplex_Stats'
            unknown: it.last().simpleName == 'Top_Unknown_Barcodes'
        }
        | set { barcodes }

    emit:
        determined   = fastq.determined
        undetermined = fastq.undetermined
        reports      = bcl_out.reports
        logs         = bcl_out.logs
        known_barcodes   = barcodes.known
        unknown_barcodes = barcodes.unknown
}

workflow {
    // Get cohorts info
    cohorts_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.bcl), file(row.sample_sheet) ] }
        | unique
        
    // Demultiplex and extract reads
    demultiplex_bcl_files(cohorts_ch)
}
