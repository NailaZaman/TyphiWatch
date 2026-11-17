/*
 * De novo assembly (SPAdes, isolate mode) and assembly statistics (QUAST)
 */

process SPADES {
    tag "$sample"
    label 'process_high'
    conda 'bioconda::spades=3.15.5'
    container params.containers.spades
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: '*.{fasta,log}'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fasta"), emit: assembly
    path "${sample}.spades.log",                emit: log

    script:
    def (r1, r2) = reads
    """
    spades.py \\
        --isolate \\
        -1 ${r1} -2 ${r2} \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        -o spades_out

    cp spades_out/contigs.fasta ${sample}.fasta
    cp spades_out/spades.log    ${sample}.spades.log
    """
}

process QUAST {
    tag "$sample"
    label 'process_low'
    conda 'bioconda::quast=5.2.0'
    container params.containers.quast
    publishDir "${params.outdir}/assemblies/quast", mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("quast_${sample}/report.tsv"), emit: report
    path "quast_${sample}",                                emit: dir

    script:
    """
    quast.py ${assembly} \\
        --output-dir quast_${sample} \\
        --labels ${sample} \\
        --min-contig 500 \\
        --threads ${task.cpus} \\
        --no-icarus --no-html --silent
    """
}
