/*
 * Read quality control: FastQC (raw), fastp (trim), MultiQC (aggregate)
 */

process FASTQC {
    tag "$sample"
    label 'process_low'
    conda 'bioconda::fastqc=0.12.1'
    container params.containers.fastqc
    publishDir "${params.outdir}/qc/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process FASTP {
    tag "$sample"
    label 'process_low'
    conda 'bioconda::fastp=0.23.4'
    container params.containers.fastp
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_R{1,2}.trimmed.fastq.gz"), emit: reads
    path "${sample}.fastp.json",                                   emit: json
    path "${sample}.fastp.html",                                   emit: html

    script:
    def (r1, r2) = reads
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${sample}_R1.trimmed.fastq.gz -O ${sample}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --cut_tail \\
        --length_required ${params.min_read_length} \\
        --thread ${task.cpus} \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html
    """
}

process MULTIQC {
    label 'process_low'
    conda 'bioconda::multiqc=1.21'
    container params.containers.multiqc
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc --force --title "TyphiWatch QC" --filename multiqc_report.html .
    """
}
