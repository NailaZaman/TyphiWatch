/*
 * Antimicrobial resistance: NCBI AMRFinderPlus on assemblies.
 * With --organism Salmonella, AMRFinderPlus also screens species-specific point
 * mutations (e.g. gyrA, parC, acrB), which covers the PointFinder-style
 * chromosomal-mutation step.
 */

process AMRFINDER_UPDATE {
    label 'process_low'
    conda 'bioconda::ncbi-amrfinderplus=3.12.8'
    container params.containers.amrfinder
    storeDir "${params.outdir}/db/amrfinderplus"      // cached: runs once

    output:
    path "amrfinder_db", emit: db

    script:
    """
    mkdir -p amrfinder_db
    amrfinder_update --force_update --database amrfinder_db
    """
}

process AMRFINDER {
    tag "$sample"
    label 'process_low'
    conda 'bioconda::ncbi-amrfinderplus=3.12.8'
    container params.containers.amrfinder
    publishDir "${params.outdir}/resistance", mode: 'copy'

    input:
    tuple val(sample), path(assembly)
    path db

    output:
    path "${sample}.amrfinder.tsv", emit: tsv

    script:
    """
    if [ -e ${db}/latest ]; then DB=${db}/latest; else DB=${db}; fi

    amrfinder \\
        --nucleotide ${assembly} \\
        --organism ${params.amr_organism} \\
        --plus \\
        --name ${sample} \\
        --database \$DB \\
        --threads ${task.cpus} \\
        --output ${sample}.amrfinder.tsv
    """
}
