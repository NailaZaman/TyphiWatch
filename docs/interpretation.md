# Interpreting TyphiWatch outputs

## Where things are
| Path | Content |
|------|---------|
| `results/qc/multiqc_report.html` | read QC, trimming and assembly statistics |
| `results/qc/failed_assembly_qc.tsv` | samples excluded by the assembly gate |
| `results/assemblies/` | contigs per sample |
| `results/typing/` | MLST, SeqSero2, Mykrobe/GenoTyphi raw outputs |
| `results/resistance/*.amrfinder.tsv` | raw AMRFinderPlus hits |
| `results/resistance/resistance_summary.tsv` | per-sample AMR calls (from `resistance_summary.py`) |
| `results/resistance/sample_summary.tsv` | AMR + typing + metadata |
| `results/phylogeny/` | core alignment, SNP distances, ML tree |
| `results/figures/` | figures and optional interactive HTML |

## Resistance categories
Rules in `scripts/resistance_summary.py` (edit there to change definitions):

| Call | Rule |
|------|------|
| Chloramphenicol | `cat*`, `cml*`, `floR`, or class PHENICOL |
| Ampicillin | any `bla*` gene |
| Co-trimoxazole | a `dfr*` gene **and** a `sul*` gene |
| Fluoroquinolone non-susceptible (FQ-NS) | any *gyrA/gyrB/parC/parE* mutation or `qnr*` gene |
| Ceftriaxone / 3GC | `blaCTX-M*`, `blaCMY*`, `blaDHA*`, or subclass CEPHALOSPORIN |
| Azithromycin | *acrB* R717 mutation, `mph(A)`, `ermB` |
| **MDR** | chloramphenicol + ampicillin + co-trimoxazole markers |
| **XDR** | MDR + FQ-NS + 3GC markers |

The XDR definition follows the usage in the Sindh, Pakistan outbreak literature
(MDR plus fluoroquinolone and third-generation-cephalosporin resistance).

`amr_profile` takes the first matching label in the order XDR > MDR > FQ-NS > Other AMR
markers > No known AMR markers.

### Cautions
* Calls are **genotype predictions**. A marker can be present without expression, and
  novel mechanisms are missed. Confirm clinically relevant results phenotypically.
* Short-read assemblies can fragment plasmids and repeat regions, so gene copy number
  and genomic context (plasmid vs chromosome) are not resolved.
* The `bla*` rule treats any beta-lactamase as ampicillin resistance; check specific
  genes in `*.amrfinder.tsv`.
* A high count of QRDR mutations (e.g. several in *gyrA/parC*) generally indicates
  stronger fluoroquinolone resistance than a single *gyrA* mutation; see
  `qrdr_mutation_count`.

## Genotypes and sequence types
* GenoTyphi genotypes come from the Mykrobe *typhi* panel. 4.3.1 is the H58 lineage.
  Check the exact nomenclature (including sublineage labels) against the current GenoTyphi
  documentation.
* The Mykrobe JSON layout changes between versions. If `genotype` is `NA` for samples
  that ran successfully, open one `*.mykrobe.json` and adjust `parse_mykrobe()`.
* MLST for *S.* Typhi is coarse (mostly ST1/ST2); genotype and SNP distance are far
  more informative.

## Phylogeny and SNP distances
* Tips are samples; the reference appears as `Reference`.
* Interpret clusters using bootstrap support (>= 95 for ultrafast bootstrap) and SNP
  distances together. There is no universal SNP threshold for "same outbreak": it
  depends on the lineage, sampling window and mutation rate, so compare against
  published thresholds for your setting.
* Long branches or isolates with many missing sites (see `core.txt`) may reflect low
  coverage or contamination, not biology.
* The tree is midpoint-rooted for display only.

## Quality flags to check first
1. Samples in `failed_assembly_qc.tsv`.
2. MultiQC: low Q30, adapter content, very low coverage.
3. `core.txt`: fraction of the reference covered per sample.
4. Species/serotype mismatches (SeqSero2 not reporting Typhi).
