#!/usr/bin/env python3
"""
TyphiWatch - figures.

Inputs : results/resistance/sample_summary.tsv, gene_matrix.tsv,
         results/phylogeny/typhiwatch.treefile (optional)
Outputs: results/figures/*.png and, if plotly is installed, an interactive HTML.

Dependencies: pandas, numpy, matplotlib (plotly optional). The Newick parser and
midpoint rooting are implemented here, so Biopython is not required.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.setrecursionlimit(20000)

PROFILE_ORDER = ["XDR", "MDR", "FQ-NS", "Other AMR markers", "No known AMR markers"]
PROFILE_COLORS = {
    "XDR": "#b2182b", "MDR": "#ef8a62", "FQ-NS": "#fddbc7",
    "Other AMR markers": "#d1e5f0", "No known AMR markers": "#67a9cf",
}


# ----------------------------------------------------------------------------
# Minimal Newick handling
# ----------------------------------------------------------------------------
class Node:
    def __init__(self):
        self.name = ""
        self.length = 0.0
        self.children = []

    @property
    def is_leaf(self):
        return not self.children


def parse_newick(text):
    text = text.strip().rstrip(";")

    def parse(i):
        node = Node()
        if text[i] == "(":
            i += 1
            while True:
                child, i = parse(i)
                node.children.append(child)
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
                raise ValueError(f"Malformed Newick near position {i}")
        j = i
        while j < len(text) and text[j] not in ",():;":
            j += 1
        node.name = text[i:j].strip("'\"")
        i = j
        if i < len(text) and text[i] == ":":
            i += 1
            j = i
            while j < len(text) and text[j] not in ",();":
                j += 1
            node.length = float(text[i:j] or 0)
            i = j
        return node, i

    root, _ = parse(0)
    return root


def midpoint_root(root):
    """Return a new rooted tree, re-rooted at the midpoint of the longest tip-to-tip path."""
    adj = {}
    names = {}
    counter = [0]

    def walk(node, parent_id):
        nid = counter[0]
        counter[0] += 1
        adj.setdefault(nid, [])
        names[nid] = node.name if node.is_leaf else ""
        if parent_id is not None:
            adj[nid].append((parent_id, node.length))
            adj[parent_id].append((nid, node.length))
        for c in node.children:
            walk(c, nid)
        return nid

    walk(root, None)
    leaves = [n for n in adj if names[n]]
    if len(leaves) < 3:
        return root

    def farthest(src):
        dist, prev, stack = {src: 0.0}, {src: None}, [src]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if v not in dist:
                    dist[v], prev[v] = dist[u] + w, u
                    stack.append(v)
        best = max(leaves, key=lambda n: dist[n])
        return best, dist, prev

    a, _, _ = farthest(leaves[0])
    b, dist_a, prev_a = farthest(a)
    half = dist_a[b] / 2.0

    # Walk from b back toward a until the cumulative distance passes the midpoint
    node = b
    while prev_a[node] is not None and dist_a[prev_a[node]] > half:
        node = prev_a[node]
    u, v = node, prev_a[node]          # edge (u, v); midpoint lies on it
    if v is None:
        return root
    w = dict(adj[u])[v]
    len_to_u = max(dist_a[u] - half, 0.0)
    len_to_v = max(w - len_to_u, 0.0)

    def build(nid, parent):
        n = Node()
        n.name = names[nid]
        for nb, wt in adj[nid]:
            if nb == parent:
                continue
            child = build(nb, nid)
            child.length = wt
            n.children.append(child)
        return n

    def build_side(start, parent, length):
        child = build(start, parent)
        child.length = length
        return child

    new_root = Node()
    new_root.children = [build_side(u, v, len_to_u), build_side(v, u, len_to_v)]
    return new_root


def ladderize(node):
    for c in node.children:
        ladderize(c)
    node.children.sort(key=lambda c: len(leaves_of(c)))
    return node


def leaves_of(node):
    if node.is_leaf:
        return [node]
    out = []
    for c in node.children:
        out.extend(leaves_of(c))
    return out


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def color_map(values, cmap_name="tab20"):
    uniq = [v for v in pd.unique(values) if pd.notna(v)]
    cmap = plt.get_cmap(cmap_name)
    return {v: cmap(i % cmap.N) for i, v in enumerate(sorted(uniq, key=str))}


def fig_profile_summary(df, out):
    counts = df["amr_profile"].value_counts().reindex(PROFILE_ORDER).fillna(0)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.barh(counts.index[::-1], counts.values[::-1], color=[PROFILE_COLORS[p] for p in counts.index[::-1]])
    for i, v in enumerate(counts.values[::-1]):
        ax.text(v + 0.05 * max(counts.max(), 1), i, int(v), va="center", fontsize=9)
    ax.set_xlabel("Isolates")
    ax.set_title("Predicted AMR profile (genotype-based)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_marker_prevalence(df, out):
    flags = {
        "Chloramphenicol": "chloramphenicol",
        "Ampicillin (bla)": "ampicillin",
        "Co-trimoxazole (dfr+sul)": "trimethoprim_sulfamethoxazole",
        "Fluoroquinolone NS (QRDR/qnr)": "fluoroquinolone_ns",
        "3rd-gen. cephalosporin (ESBL)": "ceftriaxone_3gc",
        "Azithromycin (acrB R717/mph/erm)": "azithromycin",
    }
    pct = {label: 100.0 * df[col].astype(str).str.lower().eq("true").mean() for label, col in flags.items()}
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.barh(list(pct.keys())[::-1], list(pct.values())[::-1], color="#4d4d4d")
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of isolates with resistance marker(s)")
    ax.set_title(f"Resistance marker prevalence (n={len(df)})")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_genotype_by_country(df, out):
    if "genotype" not in df.columns or df["genotype"].isna().all():
        return False
    d = df.copy()
    d["country"] = d.get("country", pd.Series("NA", index=d.index)).fillna("NA")
    top = d["genotype"].fillna("NA").value_counts().head(12).index
    d["genotype_top"] = d["genotype"].where(d["genotype"].isin(top), "other").fillna("NA")
    tab = pd.crosstab(d["country"], d["genotype_top"])
    cmap = color_map(tab.columns)
    fig, ax = plt.subplots(figsize=(8, 4))
    bottom = np.zeros(len(tab))
    for g in tab.columns:
        ax.bar(tab.index, tab[g], bottom=bottom, label=g, color=cmap[g])
        bottom += tab[g].values
    ax.set_ylabel("Isolates")
    ax.set_title("GenoTyphi genotypes by country")
    ax.legend(title="Genotype", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return True


def fig_gene_heatmap(matrix, profile, out):
    if matrix.shape[1] == 0:
        return False
    order = profile.reindex(matrix.index).map({p: i for i, p in enumerate(PROFILE_ORDER)}).fillna(99)
    matrix = matrix.loc[order.sort_values(kind="stable").index]
    keep = matrix.sum(axis=0).sort_values(ascending=False).index
    matrix = matrix[keep]
    h = max(3, 0.18 * len(matrix) + 1.5)
    w = max(5, 0.35 * matrix.shape[1] + 2.5)
    fig, ax = plt.subplots(figsize=(w, h))
    ax.imshow(matrix.values, aspect="auto", cmap="Greys", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=7)
    if len(matrix) <= 80:
        ax.set_yticks(range(len(matrix)))
        ax.set_yticklabels(matrix.index, fontsize=6)
    else:
        ax.set_yticks([])
    ax.set_title("AMR marker presence/absence")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return True


def fig_tree(tree_path, summary, out, midpoint=True):
    root = parse_newick(Path(tree_path).read_text())
    if midpoint:
        root = midpoint_root(root)
    ladderize(root)
    y, depth = {}, {}

    def assign(node, d):
        depth[id(node)] = d
        if node.is_leaf:
            y[id(node)] = len(y_leaf)
            y_leaf.append(node)
            return y[id(node)]
        ys = [assign(c, d + c.length) for c in node.children]
        y[id(node)] = sum(ys) / len(ys)
        return y[id(node)]

    y_leaf = []
    assign(root, 0.0)

    n = len(y_leaf)
    max_d = max(depth[id(t)] for t in y_leaf) or 1.0
    fig_h = max(4, 0.16 * n + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    def draw(node):
        if not node.children:
            return
        ys = [y[id(c)] for c in node.children]
        ax.plot([depth[id(node)]] * 2, [min(ys), max(ys)], color="k", lw=0.8)
        for c in node.children:
            ax.plot([depth[id(node)], depth[id(c)]], [y[id(c)]] * 2, color="k", lw=0.8)
            draw(c)

    draw(root)

    meta = summary.set_index("sample") if "sample" in summary.columns else pd.DataFrame()
    strips = [("amr_profile", "AMR profile", PROFILE_COLORS)]
    for col, label in (("genotype", "Genotype"), ("country", "Country")):
        if col in meta.columns and meta[col].notna().any():
            strips.append((col, label, color_map(meta[col])))

    strip_w = max_d * 0.04
    x0 = max_d * 1.03
    for k, (col, label, cmap) in enumerate(strips):
        for t in y_leaf:
            if t.name in meta.index:
                v = meta.at[t.name, col]
                ax.add_patch(plt.Rectangle((x0 + k * strip_w * 1.15, y[id(t)] - 0.4), strip_w, 0.8,
                                           color=cmap.get(v, "#dddddd"), lw=0))
        ax.text(x0 + k * strip_w * 1.15 + strip_w / 2, -1.2, label, rotation=90, ha="center", va="bottom", fontsize=7)

    if n <= 70:
        xl = x0 + len(strips) * strip_w * 1.15 + strip_w * 0.3
        for t in y_leaf:
            ax.text(xl, y[id(t)], t.name, va="center", fontsize=5.5)

    ax.set_ylim(n, -2.5)
    ax.set_xlim(0, x0 + len(strips) * strip_w * 1.15 + (max_d * 0.25 if n <= 70 else strip_w))
    ax.set_yticks([])
    ax.set_xlabel("Substitutions per site")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title("Core-SNP phylogeny (midpoint-rooted)" if midpoint else "Core-SNP phylogeny")

    # Legend for AMR profile
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for p, c in PROFILE_COLORS.items()]
    ax.legend(handles, PROFILE_COLORS.keys(), title="AMR profile", loc="lower left", fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def interactive_html(summary, out):
    try:
        import plotly.express as px
    except ImportError:
        print("plotly not installed - skipping interactive HTML (pip install plotly)")
        return False
    d = summary.copy()
    d["year"] = pd.to_numeric(d.get("year"), errors="coerce")
    d = d.dropna(subset=["year"])
    if d.empty:
        print("No collection years in metadata - skipping interactive HTML")
        return False
    fig = px.histogram(d, x="year", color="amr_profile", barmode="stack",
                       category_orders={"amr_profile": PROFILE_ORDER},
                       color_discrete_map=PROFILE_COLORS,
                       hover_data=[c for c in ("sample", "country", "genotype", "st") if c in d.columns],
                       title="TyphiWatch: predicted AMR profile by year")
    fig.write_html(out, include_plotlyjs="cdn")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default="results/resistance/sample_summary.tsv")
    ap.add_argument("--matrix", default="results/resistance/gene_matrix.tsv")
    ap.add_argument("--tree", default="results/phylogeny/typhiwatch.treefile")
    ap.add_argument("--outdir", default="results/figures")
    ap.add_argument("--no-midpoint", action="store_true", help="do not midpoint-root the tree")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(args.summary, sep="\t", dtype=str)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)

    fig_profile_summary(summary, out / "amr_profile_summary.png")
    fig_marker_prevalence(summary, out / "amr_marker_prevalence.png")
    fig_gene_heatmap(matrix, summary.set_index("sample")["amr_profile"], out / "amr_gene_heatmap.png")
    fig_genotype_by_country(summary, out / "genotype_by_country.png")
    interactive_html(summary, out / "typhiwatch_interactive.html")

    if Path(args.tree).exists():
        fig_tree(args.tree, summary, out / "phylogeny_annotated.png", midpoint=not args.no_midpoint)
    else:
        print(f"Tree not found ({args.tree}) - skipping phylogeny figure")

    print(f"Figures written to {out}/")


if __name__ == "__main__":
    main()
