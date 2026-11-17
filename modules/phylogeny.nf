/*
 * Comparative genomics: reference-based SNP calling (Snippy), core SNP
 * alignment (snippy-core + snp-sites), pairwise SNP distances (snp-dists)
 * and maximum-likelihood tree (IQ-TREE).
 */

process SNIPPY {
    tag "$sample"
    label 'process_medium'
    conda 'bioconda::snippy=4.6.0'
    container params.containers.snippy

    input:
    tuple val(sample), path(reads)
    path reference

    output:
    tuple val(sample), path(sample), emit: dir

    script:
    def (r1, r2) = reads
    """
    snippy \\
        --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga()} \\
        --mincov ${params.snippy_mincov} \\
        --minfrac ${params.snippy_minfrac} \\
        --ref ${reference} \\
        --R1 ${r1} --R2 ${r2} \\
        --outdir ${sample}
    """
}

process SNIPPY_CORE {
    label 'process_low'
    conda 'bioconda::snippy=4.6.0'
    container params.containers.snippy
    publishDir "${params.outdir}/phylogeny", mode: 'copy', pattern: 'core.{aln,txt,tab,vcf}'

    input:
    path sample_dirs

    output:
    path "core.aln", emit: aln
    path "core.txt", emit: stats
    path "core.tab", emit: tab

    script:
    """
    FIRST=\$(echo ${sample_dirs} | cut -d' ' -f1)
    snippy-core --ref \${FIRST}/reference/ref.fa --prefix core ${sample_dirs}
    """
}

process SNP_SITES {
    label 'process_low'
    conda 'bioconda::snp-sites=2.5.1'
    container params.containers.snp_sites
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path aln

    output:
    path "core.snps.aln", emit: aln

    script:
    """
    snp-sites -c -o core.snps.aln ${aln}
    """
}

process SNP_DISTS {
    label 'process_low'
    conda 'bioconda::snp-dists=0.8.2'
    container params.containers.snp_dists
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path aln

    output:
    path "snp_distances.tsv", emit: matrix

    script:
    """
    snp-dists -j ${task.cpus} ${aln} > snp_distances.tsv
    """
}

process IQTREE {
    label 'process_medium'
    conda 'bioconda::iqtree=2.3.4'
    container params.containers.iqtree
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path aln

    output:
    path "typhiwatch.treefile", emit: tree
    path "typhiwatch.iqtree",   emit: report
    path "typhiwatch.log",      emit: log

    script:
    """
    IQTREE=\$(command -v iqtree2 || command -v iqtree)
    NSEQ=\$(grep -c '>' ${aln})
    BOOT=""
    if [ "\$NSEQ" -ge 4 ]; then BOOT="-B ${params.iqtree_bootstraps}"; fi

    \$IQTREE -s ${aln} -m ${params.iqtree_model} \$BOOT -T ${task.cpus} --prefix typhiwatch -redo
    """
}
