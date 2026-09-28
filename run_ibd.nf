#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { GENOME }  from '../modules/plink/genome.nf'

workflow run_ibd {
    take:
    genotypes

    main:
    // Calculate pairwise IBD estimates for all samples in the merged dataset
    genotypes
        | GENOME
        | set { ibd }

    emit:
    ibd    = ibd
}
