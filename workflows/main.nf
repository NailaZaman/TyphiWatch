#!/usr/bin/env nextflow
/*
 * TyphiWatch - Salmonella Typhi WGS surveillance workflow
 *
 * Usage (from the repository root):
 *   python scripts/metadata_processing.py resolve        # build samplesheet from accessions
 *   python scripts/metadata_processing.py reference      # download CT18 reference
 *   nextflow run workflows/main.nf -profile docker
 */

nextflow.enable.dsl = 2

include { FASTQC; FASTP; MULTIQC }                                from '../modules/qc.nf'
include { SPADES; QUAST }                                         from '../modules/assembly.nf'
include { MLST; SEQSERO2; MYKROBE }                               from '../modules/typing.nf'
include { AMRFINDER_UPDATE; AMRFINDER }                           from '../modules/resistance.nf'
include { SNIPPY; SNIPPY_CORE; SNP_SITES; SNP_DISTS; IQTREE }     from '../modules/phylogeny.nf'

workflow {

    // ------------------------------------------------------------------
    // 0. Input: samplesheet.csv with columns sample,fastq_1,fastq_2
    //    (fastq_* may be local paths or http/https/ftp URLs)
    // ------------------------------------------------------------------
    ch_reads = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if( !row.sample || !row.fastq_1 || !row.fastq_2 )
                error "Samplesheet rows need sample, fastq_1 and fastq_2 (offending row: ${row})"
            tuple(row.sample, [ file(row.fastq_1), file(row.fastq_2) ])
        }

    // ------------------------------------------------------------------
    // 1. Read QC
    // ------------------------------------------------------------------
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // ------------------------------------------------------------------
    // 2. Assembly + assembly QC gate
    // ------------------------------------------------------------------
    SPADES(FASTP.out.reads)
    QUAST(SPADES.out.assembly)

    QUAST.out.report
        .map { sample, report ->
            def m = [:]
            report.readLines().each { line ->
                def p = line.split('\t')
                if( p.size() > 1 ) m[p[0]] = p[1]
            }
            def len = (m['Total length'] ?: '0') as long
            def n50 = (m['N50'] ?: '0') as long
            tuple(sample, len, n50)
        }
        .branch {
            pass: it[1] >= params.min_genome_size && it[1] <= params.max_genome_size && it[2] >= params.min_n50
            fail: true
        }
        .set { ch_gate }

    // Record samples that failed the gate
    ch_gate.fail
        .map { "${it[0]}\t${it[1]}\t${it[2]}" }
        .collectFile(
            name: 'failed_assembly_qc.tsv',
            storeDir: "${params.outdir}/qc",
            seed: "sample\ttotal_length\tn50",
            newLine: true
        )

    ch_pass = ch_gate.pass.map { tuple(it[0], true) }

    ch_assemblies = SPADES.out.assembly.join(ch_pass).map { sample, asm, flag -> tuple(sample, asm) }
    ch_trimmed    = FASTP.out.reads.join(ch_pass).map     { sample, rds, flag -> tuple(sample, rds) }

    // MultiQC over everything produced so far
    MULTIQC(
        FASTQC.out.reports.flatten()
            .mix(FASTP.out.json)
            .mix(QUAST.out.dir)
            .collect()
    )

    // ------------------------------------------------------------------
    // 3A. Genotyping
    // ------------------------------------------------------------------
    MLST(ch_assemblies)
    if( params.run_seqsero2 )  SEQSERO2(ch_assemblies)
    if( params.run_genotyphi ) MYKROBE(ch_trimmed)

    // ------------------------------------------------------------------
    // 3B. Resistome
    // ------------------------------------------------------------------
    ch_amr_db = params.amrfinder_db
        ? Channel.value(file(params.amrfinder_db, checkIfExists: true))
        : AMRFINDER_UPDATE().db
    AMRFINDER(ch_assemblies, ch_amr_db)

    // ------------------------------------------------------------------
    // 4. Comparative genomics
    // ------------------------------------------------------------------
    if( params.run_phylogeny ) {
        ch_ref = Channel.value(file(params.reference, checkIfExists: true))

        SNIPPY(ch_trimmed, ch_ref)

        ch_dirs = SNIPPY.out.dir
            .map { sample, dir -> dir }
            .collect()
            .filter { it.size() >= params.min_phylo_samples }

        SNIPPY_CORE(ch_dirs)
        SNP_SITES(SNIPPY_CORE.out.aln)
        SNP_DISTS(SNP_SITES.out.aln)
        IQTREE(SNP_SITES.out.aln)
    }
}
