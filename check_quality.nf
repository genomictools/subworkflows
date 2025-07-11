#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { FASTQC }      from '../modules/fastqc.nf'
include { SEQSTATS }    from '../modules/seqstats.nf'
include { MULTIQC }     from '../modules/multiqc.nf'

// Check quality of fastq files
workflow check_quality {
    take: 
        fastq

    main:
    // Check quality of fastq files
    fastq
        | FASTQC
        | groupTuple(by: [0, 1])
        | MULTIQC
        | set { multiqc }

    // seqkit stats
    fastq
        | SEQSTATS
        | map { it.last() }
        | collectFile (
            keepHeader: true,
            storeDir: "${params.output_dir}/seqstats",
            name: "seqstats_summary.tsv"
        )
        | set { seqstats }

    emit:
        multiqc
        seqstats
}
// [pheno, Undetermined, pheno.Undetermined.ACGATTGCTG.R1.paired, R1, /data/rds/DGE/DUDGE/MOPOPGEN/mahmed03/pipelines/demultiplex-bcl-files/test/work/9e/e366817d7c9dbc8a9ca60e54b5d4e9/pheno.Undetermined.ACGATTGCTG.R1.paired.fastq.gz]

workflow {
    // Get undetermined from cohorts
    fastq = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.sample_id, row.sample, row.read, file(row.fastq) ] }

    // Check quality of fastq files
    check_quality(fastq)
}