#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { DEEPMVP }     from '../modules/deepmvp.nf'
include { ATSNP }       from '../modules/atsnp.nf'
include { ALPHAGENOME } from '../modules/alphagenome.nf'
include { SPLICEAI }    from '../modules/spliceai.nf'
include { PANGOLIN }    from '../modules/pangolin.nf'
include { VEP }         from '../modules/vep.nf'

include { CSQ }         from '../modules/csq.nf'

include { FORMAT }      from '../modules/format.nf'
include { RESHAPE }     from '../modules/reshape.nf'
include { CONCATINATE } from '../modules/concatinate.nf'

workflow run_tool {
    take:
    variants
    tools

    main:
    variants | combine( tools ) | filter { it.last() == 'vep' }         | VEP
    variants | combine( tools ) | filter { it.last() == 'spliceai' }    | SPLICEAI
    variants | combine( tools ) | filter { it.last() == 'pangolin' }    | PANGOLIN
    variants | combine( tools ) | filter { it.last() == 'atsnp' }       | ATSNP
    variants | CSQ | filter { it.last().toInteger() > 1 }
             | combine( tools ) | filter { it.last() == 'deepmvp' }     | DEEPMVP
    variants | combine( tools ) | filter { it.last() == 'alphagenome' } | ALPHAGENOME

    // Concatenate all annotations
    DEEPMVP.out
        | concat(ALPHAGENOME.out)
        | FORMAT
        | concat(ATSNP.out)
        | RESHAPE
        // | concat(PANGOLIN.out)
        | concat(SPLICEAI.out)
        | concat(VEP.out)
        | groupTuple(by: [0,1,2])
        | CONCATINATE
        | groupTuple(by: [0,2])
        | set { annotations }

    emit:
    annotations
}

// Workflow
workflow {
    // Define input from file
    variants_ch = Channel.fromPath(params.variants)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.id, file(row.file), file(row.index) ]}
    tools_ch = Channel.from( params.tools.split(',') )

    run_tool(variants_ch, tools_ch)
}
