#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { REMOVE }      from '../modules/remove.nf'
include { FIX }         from '../modules/fix.nf'
include { FILL }        from '../modules/fill.nf'
include { CONVERT }     from '../modules/convert.nf'
include { PRUNE }       from '../modules/prune.nf'
include { COMBINE }     from '../modules/combine.nf'

fasta       = Channel.fromFilePairs(params.fasta, flat: true)
ld_regions  = Channel.fromPath(params.ld_regions)

workflow prepare_variants {
    take:
    variants
    population
    
    main:
    // Remove and fix
    variants
        | ( params.remove ? REMOVE : map {it} )
        | filter { it.last().toInteger() > 0 }
        | ( params.fix    ? combine(fasta) : map {it} )
        | ( params.fix    ? FIX    : map {it} )
        | filter { it.last().toInteger() > 0 }
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Fill and convert
    snps.cases
        | ( params.fill ? FILL : map {it} )
        | filter { it.last().toInteger() > 0 }
        | concat(snps.references)
        | CONVERT
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Prune and combine
    snps.cases
        | ( params.prune ? combine(ld_regions) : map {it} )
        | ( params.prune ? PRUNE : map {it} )
        | concat(snps.references)
        | groupTuple(by: [0,1])
        | COMBINE
        | combine(population, by: 0)
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    emit:
    cases      = snps.cases
    references = snps.references
}

// worflow
workflow {
    variants_ch = Channel.fromPath(params.variants)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ 
            row.cohort, row.type, row.chrom, row.chunk,
            file(row.vcf), filw(row.index), row.n_vars
         ] }

    population_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.population) ] }

    prepare_variants(variants_ch)
}