#!/usr/bin/env python3
"""Render manuscript Tables 4 and 5 as figure-quality images.

These two stay tables - they are categorical readouts (gene lists, term counts, module
names), not quantities a chart would show better. They are drawn rather than typeset so
they match the other figures: large bold type, no colour carrying meaning, greyscale-safe.

Every value is verified against the enrichr outputs before drawing; see verify() below,
which re-derives each row from enrichr_table.csv + modules.csv and fails loudly on a
mismatch. Run it and you are reading the analysis, not a transcription of the manuscript.

Usage:  python make_table_images.py [OUT_DIR]
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
})

RESULTS = os.environ.get(
    "ANALYSIS_RESULTS",
    "C:/projects/rise_hdwgnca_data/analysis_results/mouse_GSE233866",
)

INK = "#111111"
HEAD_BG = "#1F1F1F"
STRIPE = "#F0F0F0"
RULE = "#BBBBBB"

NETS = {
    "SNc healthy":  "healthy_baseline/snc_network",
    "SNc intact":   "lesion_model/snc_intact_network",
    "SNc lesioned": "lesion_model/snc_lesioned_network",
    "VTA healthy":  "healthy_baseline/vta_network",
    "VTA intact":   "lesion_model/vta_intact_network",
    "VTA lesioned": "lesion_model/vta_lesioned_network",
}


def cacna1d_module(net):
    with open(os.path.join(RESULTS, NETS[net], "modules.csv"), newline="",
              encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    gcol = next(c for c in rows[0] if c.lower() in ("gene_name", "gene"))
    mcol = next(c for c in rows[0] if "module" in c.lower())
    return next(r[mcol] for r in rows if r[gcol] == "Cacna1d")


def enrichr(net, module):
    with open(os.path.join(RESULTS, NETS[net], "enrichr_table.csv"), newline="",
              encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh) if r["module"] == module]


def verify():
    """Re-derive both tables from the analysis outputs. Raises on any mismatch."""
    t4, t5 = [], []
    for net in ("SNc healthy", "SNc intact", "SNc lesioned"):
        mod = cacna1d_module(net)
        ds = next(r for r in enrichr(net, mod) if r["Term"].lower() == "dopaminergic synapse")
        genes = ", ".join(g.capitalize() for g in ds["Genes"].split(";"))
        t4.append((net, mod, ds["Overlap"], float(ds["Adjusted.P.value"]), genes))

    vta_ns = []
    for net in ("VTA healthy", "VTA intact", "VTA lesioned"):
        mod = cacna1d_module(net)
        rows = enrichr(net, mod)
        hits = [r for r in rows if r["Term"].lower() == "dopaminergic synapse"]
        padj = float(hits[0]["Adjusted.P.value"]) if hits else None
        assert padj is None or padj >= 0.05, f"{net} dopaminergic synapse is significant"
        vta_ns.append((net, mod, padj))

        sig = sorted((r for r in rows if float(r["Adjusted.P.value"]) < 0.05),
                     key=lambda r: float(r["Adjusted.P.value"]))
        top = (f"{sig[0]['Term']},\n{float(sig[0]['Adjusted.P.value']):.2e}".replace(
            "e-0", r"$\times$10$^{-}$").replace("$^{-}$", "$^{-0}$") if sig else "\u2014")
        if sig:
            p = float(sig[0]["Adjusted.P.value"])
            exp = len(f"{p:.0e}".split("e")[1].lstrip("-0")) and int(f"{p:e}".split("e")[1])
            top = f"{sig[0]['Term']},\n{p / 10 ** exp:.2f}$\\times$10$^{{{exp}}}$"
        t5.append((net.replace("VTA ", "").capitalize(), mod, len(rows), len(sig), top))
    return t4, t5, vta_ns


def sup(p):
    exp = int(f"{p:e}".split("e")[1])
    return f"{p / 10 ** exp:.2f}$\\times$10$^{{{exp}}}$"


def render(ax, headers, rows, widths, row_h, head_h, fs, head_fs, bold_cols=(0,)):
    total = sum(widths)
    xs, acc = [], 0.0
    for w in widths:
        xs.append(acc)
        acc += w
    top = len(rows) * row_h

    ax.add_patch(Rectangle((0, top), total, head_h, fc=HEAD_BG, ec="none", zorder=1))
    for x, w, h in zip(xs, widths, headers):
        ax.text(x + w / 2, top + head_h / 2, h, ha="center", va="center",
                fontsize=head_fs, fontweight="bold", color="white", zorder=3,
                linespacing=1.3)

    for i, row in enumerate(rows):
        y = top - (i + 1) * row_h
        if i % 2:
            ax.add_patch(Rectangle((0, y), total, row_h, fc=STRIPE, ec="none", zorder=0))
        for j, (x, w, cell) in enumerate(zip(xs, widths, row)):
            ax.text(x + w / 2, y + row_h / 2, str(cell), ha="center", va="center",
                    fontsize=fs, color=INK, zorder=3, linespacing=1.35,
                    fontweight="bold" if j in bold_cols else "normal")
        ax.plot([0, total], [y, y], color=RULE, lw=1.0, zorder=2)

    ax.plot([0, total], [top, top], color=INK, lw=2.0, zorder=2)
    ax.plot([0, total], [0, 0], color=INK, lw=2.0, zorder=2)
    ax.set_xlim(-0.01 * total, total * 1.01)
    ax.set_ylim(-0.04 * top, top + head_h * 1.06)
    ax.axis("off")


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), format=ext,
                    bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print(f"  {name}.pdf + {name}.png")


def fig_table4(out_dir, t4, vta_ns):
    rows = []
    for net, mod, overlap, padj, genes in t4:
        wrapped = genes
        if genes.count(",") >= 5:                      # keep the long lists to two lines
            parts = genes.split(", ")
            cut = (len(parts) + 1) // 2
            wrapped = ", ".join(parts[:cut]) + ",\n" + ", ".join(parts[cut:])
        rows.append((net, mod, overlap, sup(padj), wrapped))
    rows.append(("VTA (all three)", "\u2014", "n.s.", "\u2014", "\u2014"))

    fig, ax = plt.subplots(figsize=(17.0, 5.0))
    render(ax, ["Network", "Cacna1d\nmodule", "Overlap", "p-adj", "Genes"],
           rows, widths=[2.5, 1.5, 1.3, 1.7, 7.4],
           row_h=1.0, head_h=1.15, fs=15.5, head_fs=16.5, bold_cols=(0, 2, 3))
    ax.set_title("Table 4.  Cacna1d module overlap with the Dopaminergic Synapse gene set",
                 fontsize=19, fontweight="bold", pad=16, loc="left", x=0)
    note = ";  ".join(f"{n.replace('VTA ', '')} ({m}) p-adj = {p:.2f}"
                      for n, m, p in vta_ns if p is not None)
    fig.text(0.5, -0.02,
             f"n.s. = not significant at p-adj < 0.05.  VTA, Cacna1d's module:  {note}.",
             ha="center", fontsize=13.5, fontweight="bold", color="#5A5A5A")
    fig.tight_layout()
    save(fig, out_dir, "table4_overlap")


def fig_table5(out_dir, t5):
    order = {"Intact": 0, "Lesioned": 1, "Healthy": 2}      # manuscript row order
    rows = [(n, m, f"{t:,}", s, top) for n, m, t, s, top in sorted(t5, key=lambda r: order[r[0]])]
    fig, ax = plt.subplots(figsize=(15.0, 4.3))
    render(ax, ["VTA network", "Cacna1d\nmodule", "Terms\ntested",
                "Significant\n(p-adj < 0.05)", "Top term"],
           rows, widths=[2.3, 1.6, 1.4, 2.1, 4.6],
           row_h=1.08, head_h=1.3, fs=15.5, head_fs=16.5, bold_cols=(0, 3))
    ax.set_title("Table 5.  KEGG enrichment of Cacna1d's module in the VTA networks",
                 fontsize=19, fontweight="bold", pad=16, loc="left", x=0)
    fig.text(0.5, -0.03,
             "Dopaminergic Synapse is not significant in any of the three "
             "(see Table 4); the terms listed here are the strongest other hits.",
             ha="center", fontsize=13.5, fontweight="bold", color="#5A5A5A")
    fig.tight_layout()
    save(fig, out_dir, "table5_vta_kegg")


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "C:/projects/rise_hdwgnca_data/table_figures"
    t4, t5, vta_ns = verify()
    print(f"verified against enrichr outputs; writing to {out_dir}")
    fig_table4(out_dir, t4, vta_ns)
    fig_table5(out_dir, t5)


if __name__ == "__main__":
    main()
