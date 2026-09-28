#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { GENOTYPE } from '../modules/bcftools/genotype.nf'
include { CONVERT }  from '../modules/plink/convert.nf'
include { MERGE }    from '../modules/plink/merge.nf'
include { FILTER }   from '../modules/plink/filter.nf'
include { SUBSET }   from '../modules/plink/subset.nf'

workflow align_genotypes {
    take:
    cohorts

    main:
    // Branch input by file type
    def file_types = branchCriteria { it ->
        reports : it[2] == 'report'
        vcfs    : it[2] == 'vcf'
        plinks  : it[2] == 'plink'
    }

    // Process each branch accordingly
    cohorts
        | branch(file_types)
        | set { branches }

    //  Process reports into VCFs, pass through existing VCFs, and convert to PLINK format
    branches.reports
        | GENOTYPE
        | concat(branches.vcfs)
        | CONVERT
        | set { vcfs }

    // Merge all PLINK files together, then filter and subset as needed
    Channel.empty()
        | concat(branches.plinks)
        | concat(vcfs)
        | groupTuple(by: 0)
        | MERGE
        | FILTER
        | SUBSET
        | set { genotypes }

    emit:
    genotypes
}
