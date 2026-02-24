process NOVOLOCI {
    tag "$meta.id"
    label 'process_high'

    // Placeholder for novoloci container - update when available
    container 'quay.io/biocontainers/novoloci:latest'

    input:
    tuple val(meta), path(reads1), path(reads2), path(seed)

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"), emit: assembly
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "NOVOloci is not yet implemented. Placeholder for future implementation."
    touch "${meta.id}_assembly.fasta"
    touch "${meta.id}.log"
    
    cat > versions.yml << EOF
    "${task.process}":
        novoloci: "placeholder"
    EOF
    """

    stub:
    """
    touch "${meta.id}_assembly.fasta"
    touch "${meta.id}.log"
    
    cat > versions.yml << EOF
    "${task.process}":
        novoloci: "placeholder"
    EOF
    """
}
