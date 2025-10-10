#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PICK }        from '../modules/bcftools/pick.nf'
include { REMOVE }      from '../modules/bcftools/remove.nf'
include { FIX }         from '../modules/bcftools/fix.nf'
include { FILL }        from '../modules/bcftools/fill.nf'
include { CONVERT }     from '../modules/plink/convert.nf'
include { PRUNE }       from '../modules/plink/prune.nf'
include { COMBINE }     from '../modules/plink/combine.nf'

fasta       = Channel.fromFilePairs(params.fasta, flat: true)
ld_regions  = Channel.fromPath(params.ld_regions)

workflow prepare_variants {
    take:
    cohorts
    variants
    population
    
    main:
    // Remove and fix
    cohorts
        | combine(variants)
        | combine(population, by: 0)
        | PICK
        | filter { it.last().toInteger() > 0 }
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
    cohorts_ch = Channel.fromPath(params.variants)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ 
            row.cohort, row.type, row.chrom, row.chunk,
            file(row.vcf), file(row.index), row.n_vars
         ] }

    variants_ch   = Channel.fromPath(params.snplist)

    population_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, file(row.population) ] }

    prepare_variants( cohorts_ch, variants_ch, population_ch)
}
