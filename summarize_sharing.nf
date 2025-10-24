#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SHARING }     from '../modules/rvs/sharing.nf'
include { CLASSIFY }    from '../modules/rvs/classify.nf'
include { ATTACH }      from '../modules/rvs/attach.nf'
include { DRAW }        from '../modules/rvs/draw.nf'

type_ch = Channel.of( 'variant', 'gene' )

workflow summarize_sharing {
    take:
    variants
    family
    blacklist
    to_draw
    
    main:
    // Extract variants stats
    variants
        | combine(family, by: 0)
        | combine(blacklist)
        | SHARING
        | groupTuple(by: [0, 1])
        | combine(type_ch)
        | CLASSIFY
        | set { shared }

    // Draw pedigrees
    if ( params.draw ) {
        variants
            | combine(family, by: 0)
            | ATTACH
            | combine(to_draw, by: 0)
            | DRAW
    }

    emit:
    shared
}

workflow  {
    family_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.famid, file(row.cases), file(row.pedigree)] }
        | unique

    variants_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.famid, row.category, file(row.rlist), file(row.annotation)
        ] }

    blacklist_ch = Channel.fromPath(params.blacklist)

    if ( params.draw ) {
    to_draw = Channel.fromPath(params.draw_genes)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.famid, row.gene, row.variant ] }
        | groupTuple(by: [0, 1])
    } else {
    to_draw = Channel.empty()
    }

    summarize_genes( variants_ch, family_ch, blacklist_ch, to_draw )
}
