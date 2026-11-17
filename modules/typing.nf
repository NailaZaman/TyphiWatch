/*
 * Genotyping: MLST (assembly), SeqSero2 serotype (assembly), GenoTyphi genotype
 * via the Mykrobe Typhi panel (reads).
 *
 * SeqSero2 and Mykrobe use errorStrategy 'ignore' so that a failure in one
 * optional typing tool never stops the rest of the pipeline.
 */

process MLST {
    tag "$sample"
    label 'process_low'
    conda 'bioconda::mlst=2.23.0'
    container params.containers.mlst
    publishDir "${params.outdir}/typing/mlst", mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    path "${sample}.mlst.tsv", emit: tsv

    script:
    def scheme = params.mlst_scheme ? "--scheme ${params.mlst_scheme}" : ''
    """
    mlst ${scheme} --label ${sample} --threads ${task.cpus} ${assembly} > ${sample}.mlst.tsv
    """
}

process SEQSERO2 {
    tag "$sample"
    label 'process_low'
    errorStrategy 'ignore'
    conda 'bioconda::seqsero2=1.2.1'
    container params.containers.seqsero2
    publishDir "${params.outdir}/typing/seqsero2", mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    path "${sample}.seqsero2.tsv", emit: tsv

    script:
    // -t 4 = genome assembly input; -m a = allele-based mode (required for assemblies)
    """
    SeqSero2_package.py -t 4 -m a -p ${task.cpus} -i ${assembly} -n ${sample} -d seqsero2_out
    cp seqsero2_out/SeqSero_result.tsv ${sample}.seqsero2.tsv
    """
}

process MYKROBE {
    tag "$sample"
    label 'process_medium'
    errorStrategy 'ignore'
    conda 'bioconda::mykrobe=0.13.0'
    container params.containers.mykrobe
    publishDir "${params.outdir}/typing/mykrobe", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path "${sample}.mykrobe.json", emit: json

    script:
    def (r1, r2) = reads
    // NOTE: check `mykrobe predict --help` for your installed version if flags differ.
    """
    mykrobe predict \\
        --sample ${sample} \\
        --species typhi \\
        --format json \\
        --out ${sample}.mykrobe.json \\
        --threads ${task.cpus} \\
        -1 ${r1} ${r2}
    """
}
