#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { SPLITFASTA }  from '../modules/seqkit/splitfasta.nf'
include { INDEXFASTA }  from '../modules/samtools/indexfasta.nf'
include { DICTFASTA }   from '../modules/gatk/dictfasta.nf'

workflow split_fasta {
    take:
    assembly
    fasta

    main:
    Channel.fromPath(fasta) 
        | map { [ assembly, it ] }
        | SPLITFASTA
        | transpose
        | map { [it.first(), it.last().name.split('\\.')[1], it.last()] }
        | ( params.remove_nonstandard_chroms ? filter { !(it.last().name ==~ /.*(?:decoy|alt|HLA|chrUn|random).*/) } : map { it } )
        | take(3)
        | INDEXFASTA
        | DICTFASTA 
        | set { fasta_ch }

    emit:
    fasta_ch
}

workflow {
    split_fasta(params.assembly, params.fasta) 
}
