#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BCFTOOLS_SUBSET } from '../modules/subset.nf'
include { SAMTOOLS_SUBSET } from '../modules/subset.nf'
include { EXTRACT }       from '../modules/extract.nf'
include { ONCOEXTRACT }   from '../modules/oncoextract.nf'
include { STATS }         from '../modules/stats.nf'
include { PLOT }          from '../modules/plot.nf'

workflow infer_loss {
    take:
    vcf
    bam

    main:
    vcfs = vcf | BCFTOOLS_SUBSET
    bams = bam | SAMTOOLS_SUBSET

    vcf | STATS 
    vcfs
        | combine( bams, by: [0,1,2,3] )
        | ( params.tumour_only ? filter { it[3] == 'tumor' } : groupTuple(by: [0,1,2]) )
        | combine(Channel.fromPath(params.fasta))
        | ( params.tumour_only ? EXTRACT : ONCOEXTRACT )
}

// workflow
workflow {
    vcf_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ 
            row.cohort, row.key, row.sample, row.sample_type,
            file(row.vcf), file(row.vcf_index)
        ] }

    bam_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ 
            row.cohort, row.key, row.sample, row.sample_type,
            file(row.bam), file(row.bam_index)
        ] }

    infer_loss(vcf_ch, bam_ch)
}