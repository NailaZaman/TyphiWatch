# Methodology

TyphiWatch takes paired-end Illumina reads of *Salmonella enterica* serovar Typhi
isolates and produces QC metrics, assemblies, genotypes, resistance predictions and a
core-SNP phylogeny. It is implemented in Nextflow (DSL2); tools run in containers.

## 0. Input preparation (`scripts/metadata_processing.py`)
* `search` queries the ENA Portal API for Illumina paired-end runs under taxon 90370
  (*S.* Typhi), optionally filtered by country.
* `resolve` expands run / experiment / study accessions in `data/accession_list.tsv`
  into run-level records and writes `data/metadata/metadata.tsv` (country, year,
  instrument, ...) and `data/metadata/samplesheet.csv` (sample, fastq_1, fastq_2 as
  ENA https URLs; Nextflow downloads them).
* `reference` downloads the reference genome (default CT18, AL513382.1) as GenBank.

Only paired-end runs are kept. Sample IDs are the run accession (non-alphanumerics
replaced by `_`).

## 1. Read QC (`modules/qc.nf`)
* **FastQC** on raw reads.
* **fastp**: adapter removal (paired-end detection), tail trimming at Q20, minimum read
  length `params.min_read_length` (default 50).
* **MultiQC** aggregates FastQC, fastp and QUAST.

## 2. Assembly (`modules/assembly.nf`)
* **SPAdes** `--isolate` (recommended for high-coverage single-isolate data).
* **QUAST** with a 500 bp minimum contig.
* **QC gate** (in `workflows/main.nf`): samples are excluded from typing, AMR and
  phylogeny if total assembly length is outside `[min_genome_size, max_genome_size]`
  (default 4.3-5.5 Mb) or N50 < `min_n50` (default 20 kb). Excluded samples are listed in
  `results/qc/failed_assembly_qc.tsv`. A wrong genome size often indicates
  contamination, a mixed culture, low coverage or a non-Typhi isolate.

## 3A. Typing (`modules/typing.nf`)
* **MLST** (auto-detected scheme; *S.* Typhi is typically ST1 or ST2).
* **SeqSero2** (allele mode on assemblies) for serotype prediction.
* **GenoTyphi genotype** from the **Mykrobe** *typhi* panel on trimmed reads (e.g.
  4.3.1 = H58). SeqSero2 and Mykrobe are non-blocking (`errorStrategy 'ignore'`).

## 3B. Resistome (`modules/resistance.nf`)
* **AMRFinderPlus** on assemblies with `--organism Salmonella --plus`. The organism flag
  enables species-specific point-mutation screening (e.g. *gyrA*, *parC*, *acrB*), which
  covers the chromosomal-mutation role of PointFinder. A separate ResFinder/PointFinder
  run is not included.
* `scripts/resistance_summary.py` converts hits into per-sample calls (see
  `docs/interpretation.md` for the exact rules).

## 4. Comparative genomics (`modules/phylogeny.nf`)
* **Snippy**: map trimmed reads to the reference and call SNPs
  (`--mincov 10 --minfrac 0.9`).
* **snippy-core**: core alignment across isolates; **snp-sites -c** keeps variable,
  unambiguous columns.
* **snp-dists**: pairwise SNP distance matrix.
* **IQ-TREE 2**: `GTR+G+ASC` (ascertainment-bias correction, appropriate because
  constant sites were removed) with 1000 ultrafast bootstraps (skipped if < 4 sequences).
* Requires at least `min_phylo_samples` (default 4) isolates passing QC.

## 5. Reporting (`scripts/resistance_summary.py`, `scripts/visualization.py`)
Run after the workflow. Produces summary tables, an AMR marker heatmap, resistance
prevalence, genotype-by-country bars, an annotated midpoint-rooted tree and (with
plotly) an interactive HTML.

## Known limitations
* Recombination is not masked (Typhi has low recombination but mobile elements and
  phage regions can inflate SNP counts). For publication-grade trees, mask these
  regions and/or run Gubbins.
* The reference (CT18) is a single genome; mapping to it ignores accessory content.
* Rooting is by midpoint; use an outgroup for formal inference.
* AMR output is genotype-based; it is not phenotypic susceptibility.
* Public metadata (country, year) is often incomplete or inconsistently formatted.
