#!/usr/bin/env bash
# TyphiWatch preflight: checks that your machine is ready before the first run.
# Usage (from the repo root, with the conda env active):  bash scripts/preflight.sh
# Prints OK / WARN / MISSING with a hint for each check. Changes nothing.

pass=0; warn=0; fail=0
ok()   { echo "  [OK]      $1"; pass=$((pass+1)); }
warn() { echo "  [WARN]    $1"; warn=$((warn+1)); }
miss() { echo "  [MISSING] $1"; fail=$((fail+1)); }

echo "== Environment =="
if grep -qi microsoft /proc/version 2>/dev/null; then ok "running inside WSL"; else warn "not WSL (fine on native Linux/macOS)"; fi
case "$PWD" in
  /mnt/*) warn "project is on the Windows drive ($PWD); move it into ~ (Linux filesystem) for speed" ;;
  *)      ok "project is on the Linux filesystem ($PWD)" ;;
esac

mem_gb=$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo 2>/dev/null || echo 0)
if   [ "$mem_gb" -ge 20 ]; then ok "RAM visible to Linux: ${mem_gb} GB"
elif [ "$mem_gb" -ge 12 ]; then warn "RAM visible to Linux: ${mem_gb} GB (raise memory= in C:\\Users\\<you>\\.wslconfig, then 'wsl --shutdown')"
else miss "RAM visible to Linux: ${mem_gb} GB (too low for SPAdes)"; fi

free_gb=$(df -Pk . | awk 'NR==2 {printf "%d", $4/1024/1024}')
if [ "$free_gb" -ge 100 ]; then ok "free disk here: ${free_gb} GB"
elif [ "$free_gb" -ge 40 ]; then warn "free disk here: ${free_gb} GB (enough for a small run only)"
else miss "free disk here: ${free_gb} GB"; fi

echo "== Tools =="
for t in git curl unzip python3; do
  if command -v "$t" >/dev/null 2>&1; then ok "$t"; else miss "$t  -> sudo apt install -y $t"; fi
done
if command -v java >/dev/null 2>&1; then ok "java: $(java -version 2>&1 | head -1)"; else miss "java  -> activate the env: conda activate typhiwatch"; fi
if command -v nextflow >/dev/null 2>&1; then ok "nextflow: $(nextflow -version 2>&1 | grep -m1 -i version)"; else miss "nextflow  -> conda activate typhiwatch"; fi
if command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1; then ok "conda/mamba"; else warn "conda/mamba not found (needed to create the env)"; fi

echo "== Docker =="
if command -v docker >/dev/null 2>&1; then
  if docker info >/dev/null 2>&1; then ok "docker daemon reachable"
  else miss "docker installed but not reachable -> start Docker Desktop and enable WSL integration for this distro"; fi
else
  miss "docker not found in WSL -> install Docker Desktop and enable WSL integration (or use -profile conda)"
fi

echo "== Python packages =="
for m in pandas numpy matplotlib; do
  if python3 -c "import $m" 2>/dev/null; then ok "python: $m"; else miss "python: $m  -> conda activate typhiwatch"; fi
done
if python3 -c "import plotly" 2>/dev/null; then ok "python: plotly (interactive HTML)"; else warn "python: plotly missing (optional)"; fi

echo "== Network =="
if curl -fs --max-time 25 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRX7117257&result=read_run&fields=run_accession&format=tsv" 2>/dev/null | grep -q "run_accession"; then
  ok "ENA API reachable"; else miss "cannot reach the ENA API (check internet/VPN/proxy)"; fi
if curl -fs --max-time 25 -o /dev/null "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"; then ok "NCBI E-utilities reachable"; else warn "cannot reach NCBI (needed for the reference download)"; fi
if curl -fs --max-time 25 -o /dev/null -I "https://quay.io/v2/"; then ok "quay.io reachable (container images)"; else warn "cannot reach quay.io"; fi

echo "== Project files =="
for f in nextflow.config workflows/main.nf modules/qc.nf modules/assembly.nf modules/typing.nf modules/resistance.nf modules/phylogeny.nf \
         scripts/metadata_processing.py scripts/resistance_summary.py scripts/visualization.py data/accession_list.tsv; do
  if [ -f "$f" ]; then ok "$f"; else miss "$f (are you in the TyphiWatch folder?)"; fi
done

echo
echo "Summary: $pass OK, $warn warnings, $fail missing"
[ "$fail" -eq 0 ] && echo "Ready: next run  python scripts/metadata_processing.py resolve --max-per-accession 4" || echo "Fix the MISSING items first, then rerun this script."
exit 0
