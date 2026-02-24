/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NOVOPLASTY           } from '../../modules/local/novoplasty/main'
include { NOVOLOCI             } from '../../modules/local/novoloci/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY {

    take:
    ch_reads // channel: [ val(meta), val(read_type), path(read1), path(read2), path(seed) ]

    main:
    ch_versions = channel.empty()

    // Split reads by type
    ch_illumina = ch_reads.filter { it[1] == 'illumina' }
        .map { meta, read_type, read1, read2, seed -> [ meta, read1, read2, seed ] }
    
    ch_long_reads = ch_reads.filter { it[1] in ['ONT', 'Pacbio'] }
        .map { meta, read_type, read1, read2, seed -> [ meta, read1, read2, seed ] }
    
    ch_genomes = ch_reads.filter { it[1] == 'genome' }
        .map { meta, read_type, read1, read2, seed -> [ meta, read1 ] }

    // Run NOVOPlasty for illumina reads
    NOVOPLASTY(ch_illumina)
    ch_versions = ch_versions.mix(NOVOPLASTY.out.versions)

    // Run NOVOloci for long reads
    NOVOLOCI(ch_long_reads)
    ch_versions = ch_versions.mix(NOVOLOCI.out.versions)

    // Combine results from all sources
    ch_assemblies = NOVOPLASTY.out.assembly
        .concat(NOVOLOCI.out.assembly)
        .concat(ch_genomes)

    emit:
    assemblies = ch_assemblies
    versions   = ch_versions
}
