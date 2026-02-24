/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CMSEARCH_RRNA         } from '../../modules/local/cmsearch_rrna/main'
include { EXTRACT_RRNA          } from '../../modules/local/extract_rrna/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RRNA DETECTION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RRNA_DETECTION {

    take:
    ch_assemblies

    main:
    ch_versions = channel.empty()

    CMSEARCH_RRNA(ch_assemblies)
    ch_versions = ch_versions.mix(CMSEARCH_RRNA.out.versions)

    ch_for_extraction = ch_assemblies
        .join(CMSEARCH_RRNA.out.ssu_tab)
        .join(CMSEARCH_RRNA.out.lsu_tab)

    EXTRACT_RRNA(ch_for_extraction)
    ch_versions = ch_versions.mix(EXTRACT_RRNA.out.versions)

    emit:
    ssu      = EXTRACT_RRNA.out.ssu_fasta
    its      = EXTRACT_RRNA.out.its_fasta
    lsu      = EXTRACT_RRNA.out.lsu_fasta
    rrna_log      = EXTRACT_RRNA.out.rrna_log
    ssu_tab  = CMSEARCH_RRNA.out.ssu_tab
    s58_tab  = CMSEARCH_RRNA.out.s58_tab
    lsu_tab  = CMSEARCH_RRNA.out.lsu_tab
    versions = ch_versions
}