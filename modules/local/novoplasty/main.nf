process NOVOPLASTY {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/novoplasty:4.3.5--pl5321hdfd78af_0'

    input:
    tuple val(meta), path(reads1), path(reads2), path(seed)

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"), emit: assembly
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def forward_reads = reads1 ? reads1.toString() : ''
    def reverse_reads = reads2 ? reads2.toString() : ''
    def seed_input = seed.toString()
    def genome_range = params.genome_range ?: '6000-20000'
    def kmer = params.kmer ?: 33
    def max_memory = (task.memory.toGiga() as Integer) - 1
    
    """
    # Create config file from template
    cat > novoplasty_${meta.id}.config << 'CONFIGEOF'
Project:
-----------------------
Project name          = ${meta.id}
Type                  = mito
Genome Range          = ${genome_range}
K-mer                 = ${kmer}
Max memory            = ${max_memory}
Extended log          = 0
Save assembled reads  = no
Seed Input            = ${seed_input}
Extend seed directly  = no
Reference sequence    = 
Variance detection    = no
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = ${forward_reads}
Reverse reads         = ${reverse_reads}
Store Hash            = 

Heteroplasmy:
-----------------------
MAF                   = 
HP exclude list       = 
PCR-free              = 

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Reduce ambigious N's  = 
Output path           = ./novoplasty_output/
CONFIGEOF

    NOVOPlasty4.3.5.pl -c novoplasty_${meta.id}.config

    # Rename output assembly
    if [ -f "novoplasty_output/Circularized_assembly/${meta.id}_Circularized_*.fasta" ]; then
        cp "novoplasty_output/Circularized_assembly/${meta.id}_Circularized_"*.fasta "${meta.id}_assembly.fasta"
    elif [ -f "novoplasty_output/Assembled_sequences/${meta.id}_*.fasta" ]; then
        cp "novoplasty_output/Assembled_sequences/${meta.id}_"*.fasta "${meta.id}_assembly.fasta"
    else
        # Create empty file if assembly failed
        touch "${meta.id}_assembly.fasta"
    fi

    cat > versions.yml << EOF
    "${task.process}":
        novoplasty: \$(NOVOPlasty4.3.5.pl -v 2>&1 | grep -oP '(?<=NOVOPlasty)[0-9.]+' || echo "4.3.5")
    EOF
    """

    stub:
    """
    touch "${meta.id}_assembly.fasta"
    touch "${meta.id}.log"
    
    cat > versions.yml << EOF
    "${task.process}":
        novoplasty: 4.3.5
    EOF
    """
}
