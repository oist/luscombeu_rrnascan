process CMSEARCH_RRNA {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/infernal:1.1.5--pl5321h7b50bb2_4'

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${meta.id}_SSU.tab"), emit: ssu_tab
    tuple val(meta), path("${meta.id}_5.8S.tab"), emit: s58_tab
    tuple val(meta), path("${meta.id}_LSU.tab"), emit: lsu_tab
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def models_dir = "${projectDir}/assets/models"
    def args = task.ext.args ?: '-E 1e-15'
    
    """
    cmsearch $args --tblout ${meta.id}_SSU.tab  ${models_dir}/SSU_eukarya_RF01960 ${sequence} > /dev/null
    cmsearch $args --tblout ${meta.id}_5.8S.tab ${models_dir}/5_8S_rRNA_RF00002   ${sequence} > /dev/null
    cmsearch $args --tblout ${meta.id}_LSU.tab  ${models_dir}/LSU_eukarya_RF02543 ${sequence} > /dev/null

    cat > versions.yml << EOF
    "${task.process}":
        infernal: \$(cmsearch --version | grep "^cmsearch" | awk '{print \$3}')
    EOF
    """

    stub:
    """
    touch ${meta.id}_SSU.tab
    touch ${meta.id}_5.8S.tab
    touch ${meta.id}_LSU.tab
    
    cat > versions.yml << EOF
    "${task.process}":
        infernal: "1.1.5"
    EOF
    """
}
