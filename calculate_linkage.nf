#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { RECODE }    from '../modules/plink/recode.nf'
include { STATS }     from '../modules/merlin/stats.nf'
include { LINKAGE }   from '../modules/merlin/linkage.nf'

workflow calculate_linkage {
    take:
    variants
    test_ch

    main:
    // Extract variants stats
    variants    
        | RECODE
        | STATS
        | combine(test_ch)
        | LINKAGE
        | set { linkage }

    emit:
    linkage
}

workflow  {
    variants_ch = Channel.fromPath(params.cohorts)
        | splitCsv(header: true, sep: ',')
        | map { row -> [
            row.famid, row.phenotype, row.chrom, file(row.freq), file(row.bim), file(row.bed), file(row.fam)
        ] }
    
    test_ch = Channel.of(params.test.split(','))

    linkage = calculate_linkage( variants_ch, test_ch )

    linkage
        | collectFile (
            keepHeader: true,
            storeDir: "${params.output_dir}/summary",
        )
        { it -> [ "${it[0]}.${it[2]}.tsv", it.last() ] }
}
