#!/usr/bin/env python3
"""Convert manuscript Tables 6-8 into figures.

  Table 6 -> paired dot plot    kME of four prior Cav1.3 partner genes, intact vs lesioned,
                                against Cacna1d's own kME as the reference
  Table 7 -> heatmap            13 KEGG dopaminergic-synapse genes x 6 networks
  Table 8 -> quadrant scatter   worst SNc percentile vs best VTA percentile, with cutoffs

Tables 6 and 7 are read from the analysis CSVs, not transcribed. Table 8 has no source file
in analysis_results - its percentiles and rejection reasons are manuscript values and are
embedded below; change them here if the manuscript changes.

Styling follows the same eNeuro rules as make_figures.py: no top/right spines, no red/green.
Everything is additionally legible in greyscale - the heatmap prints its values in-cell, and
the other two separate their series by marker shape (filled vs open) rather than by colour.

Usage:  python make_table_figures.py [OUT_DIR]
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "axes.linewidth": 1.4,
})

RESULTS = os.environ.get(
    "ANALYSIS_RESULTS",
    "C:/projects/rise_hdwgnca_data/analysis_results/mouse_GSE233866",
)
KEGG = os.path.join(RESULTS, "kegg_dopaminergic_synapse")

INK = "#111111"
MUTED = "#6B6B6B"

# Table 6: the four prior-literature Cav1.3 partner genes, ordered by SNc intact kME.
T6_GENES = ["Kcnj6", "Drd2", "Kcnn3", "Nr4a2"]
T6_ALIAS = {"Kcnj6": "GIRK2", "Drd2": "", "Kcnn3": "SK3", "Nr4a2": "Nurr1"}

# Table 7: detected KEGG dopaminergic-synapse genes, in the manuscript's row order.
T7_GENES = ["Kcnj6", "Drd2", "Slc6a3", "Prkca", "Gria3", "Ddc", "Plcb4",
            "Slc18a2", "Itpr1", "Ppp2r2b", "Cacna1d", "Th", "Cacna1c"]
T7_NETS = ["SNc healthy", "SNc intact", "SNc lesioned",
           "VTA healthy", "VTA intact", "VTA lesioned"]

# Table 8: manuscript values - no CSV source in analysis_results.
# Trailing fields place the label: (dy in percentile units, horizontal alignment, x multiplier).
T8 = [
    ("Kcnj6",    0.1, 98.1, None,                             -9.5, "center", 1.00),
    ("Kcnd3",    0.6, 98.0, None,                             -9.5, "right",  0.93),
    ("Prkca",    0.7, 99.4, None,                            -19.0, "left",   1.08),
    ("Ddc",      1.2, 94.0, None,                             -9.5, "center", 1.00),
    ("Drd2",     0.2,  7.3, "not inverted",                    4.5, "center", 1.00),
    ("Kcnn3",    0.3,  1.6, "central in VTA\nhealthy/intact",  4.5, "center", 1.02),
    ("Slc6a3",   0.3, 28.0, "inverts in 2 of\n3 VTA networks", 4.5, "center", 1.00),
    ("Gria3",    0.7,  1.8, "central in VTA\nlesioned",        4.5, "center", 1.00),
    ("Nr4a2",    1.1, 27.4, "no inversion",                    4.5, "center", 1.00),
    ("Plcb4",    3.6,  4.8, "weaker in SNc,\nno inversion",    4.5, "center", 1.00),
    ("Slc18a2",  4.0, 27.0, "weaker in SNc",                   4.5, "center", 1.00),
    ("Itpr1",    1.9,  1.1, "central in\nboth regions",        4.5, "right",  0.88),
    ("Ppp2r2b",  2.1,  2.9, "central in\nboth regions",       16.0, "left",   1.10),
]
# No cutoff was specified for the VTA axis. Every kept gene sits at >=94 and every rejected
# one at <=28, so the line is drawn in that gap; it separates the data, it is not a threshold
# taken from the manuscript.
VTA_INVERSION_CUTOFF = 90.0


def bare(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def module_table(cond):
    path = os.path.join(KEGG, f"snc_{cond}_cacna1d_module_genes.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["gene_name"]: (float(r["kME"]), int(r["rank"])) for r in csv.DictReader(fh)}


def kegg_table():
    path = os.path.join(KEGG, "kegg_dopaminergic_synapse_snc_vs_vta.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["gene"]: r for r in csv.DictReader(fh)}


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), format=ext,
                    bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print(f"  {name}.pdf + {name}.png")


def fig_table6(out_dir):
    intact, lesioned = module_table("intact"), module_table("lesioned")
    fig, ax = plt.subplots(figsize=(12.5, 7.4))

    ref_i, ref_l = intact["Cacna1d"][0], lesioned["Cacna1d"][0]
    ax.axvspan(ref_i, ref_l, color="#C9C9C9", alpha=0.45, zorder=0)
    for x, ls in ((ref_i, "--"), (ref_l, ":")):
        ax.axvline(x, color=INK, ls=ls, lw=2.2, zorder=1)

    ys = list(range(len(T6_GENES)))[::-1]
    for y, g in zip(ys, T6_GENES):
        (ki, ri), (kl, rl) = intact[g], lesioned[g]
        ax.plot([ki, kl], [y, y], color=INK, lw=2.4, zorder=2, solid_capstyle="round")
        ax.scatter([ki], [y], s=430, marker="o", c=INK, edgecolors=INK, zorder=4)
        ax.scatter([kl], [y], s=430, marker="s", c="white", edgecolors=INK,
                   linewidths=2.6, zorder=4)
        lo, hi = (ki, kl) if ki <= kl else (kl, ki)
        ax.text(lo - 0.022, y, f"#{ri if ki <= kl else rl}", ha="right", va="center",
                fontsize=13, fontweight="bold", color=MUTED)
        ax.text(hi + 0.022, y, f"#{rl if ki <= kl else ri}", ha="left", va="center",
                fontsize=13, fontweight="bold", color=MUTED)

    ax.set_yticks(ys)
    ax.set_yticklabels([f"{g}\n({T6_ALIAS[g]})" if T6_ALIAS[g] else g for g in T6_GENES],
                       fontsize=19, fontweight="bold")
    ax.set_ylim(-0.75, len(T6_GENES) - 0.25)
    ax.set_xlim(0.26, 0.90)
    ax.set_xlabel("Module membership (kME) in Cacna1d's SNc module",
                  fontsize=17, fontweight="bold", labelpad=12)
    ax.tick_params(axis="x", labelsize=15)
    for lab in ax.get_xticklabels():
        lab.set_fontweight("bold")
    ax.grid(axis="x", color="#E4E4E4", lw=1.0, zorder=0)
    ax.set_axisbelow(True)
    bare(ax)

    ax.text(ref_l + 0.008, len(T6_GENES) - 0.42,
            f"Cacna1d\n{ref_i:.3f} / {ref_l:.3f}", fontsize=14, fontweight="bold",
            color=INK, ha="left", va="top")
    ax.set_title("All four prior Ca$_V$1.3 partner genes are more central to the module "
                 "than Cacna1d itself",
                 fontsize=19, fontweight="bold", pad=20)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=INK, mec=INK, ms=17, label="Intact"),
        Line2D([], [], marker="s", ls="", mfc="white", mec=INK, mew=2.6, ms=16, label="Lesioned"),
        Line2D([], [], color=INK, ls="--", lw=2.2, label="Cacna1d kME (intact)"),
        Line2D([], [], color=INK, ls=":", lw=2.2, label="Cacna1d kME (lesioned)"),
    ], loc="lower right", fontsize=14, frameon=True, prop={"size": 14, "weight": "bold"})
    fig.text(0.5, -0.015, "#n = rank within the module (1,640 genes intact; 855 lesioned)",
             ha="center", fontsize=13, color=MUTED, fontweight="bold")
    fig.tight_layout()
    save(fig, out_dir, "table6_dotplot")


def fig_table7(out_dir):
    rows = kegg_table()
    vals, stars = [], []
    for g in T7_GENES:
        r = rows[g]
        vals.append([float(r[f"{n} kME"]) for n in T7_NETS])
        stars.append([r[f"{n} in_module"].strip().upper() in ("TRUE", "T", "1", "YES")
                      for n in T7_NETS])

    lim = max(abs(v) for row in vals for v in row)
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    cmap = plt.get_cmap("PuOr_r")      # purple = negative, white = 0, orange = positive

    fig, ax = plt.subplots(figsize=(14.5, 11.6))
    im = ax.imshow(vals, cmap=cmap, norm=norm, aspect="auto")

    for i, g in enumerate(T7_GENES):
        for j in range(len(T7_NETS)):
            v = vals[i][j]
            r, gg, b, _ = cmap(norm(v))
            lum = 0.299 * r + 0.587 * gg + 0.114 * b
            ax.text(j, i, f"{v:+.3f}{'*' if stars[i][j] else ''}",
                    ha="center", va="center", fontsize=14.5, fontweight="bold",
                    color="white" if lum < 0.5 else INK)

    ax.set_xticks(range(len(T7_NETS)))
    ax.set_xticklabels([n.replace(" ", "\n") for n in T7_NETS],
                       fontsize=16, fontweight="bold")
    ax.set_yticks(range(len(T7_GENES)))
    ax.set_yticklabels(T7_GENES, fontsize=17, fontweight="bold")
    for lab in ax.get_yticklabels():
        if lab.get_text() == "Cacna1d":
            lab.set_color("#B3006B")

    ax.axvline(2.5, color=INK, lw=4.0)                      # SNc | VTA divide
    for k in range(len(T7_GENES) + 1):
        ax.axhline(k - 0.5, color="white", lw=1.6)
    for k in range(len(T7_NETS) + 1):
        ax.axvline(k - 0.5, color="white", lw=1.6)
    ax.axhline(T7_GENES.index("Cacna1d") - 0.5, color=INK, lw=2.6, ls=":")
    ax.axhline(T7_GENES.index("Cacna1d") + 0.5, color=INK, lw=2.6, ls=":")
    ax.set_xlim(-0.5, len(T7_NETS) - 0.5)
    ax.set_ylim(len(T7_GENES) - 0.5, -0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cb.set_label("kME toward Cacna1d's module", fontsize=16, fontweight="bold", labelpad=14)
    cb.ax.tick_params(labelsize=14)
    for lab in cb.ax.get_yticklabels():
        lab.set_fontweight("bold")

    ax.set_title("Dopaminergic-synapse genes track Cacna1d's module in SNc and invert in VTA",
                 fontsize=20, fontweight="bold", pad=24)
    fig.text(0.5, 0.012,
             "* gene is assigned to Cacna1d's module in that network.  "
             "Values are kME regardless of assignment.",
             ha="center", fontsize=14, fontweight="bold", color=MUTED)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    save(fig, out_dir, "table7_heatmap")


def fig_table8(out_dir):
    fig, ax = plt.subplots(figsize=(16.0, 9.8))

    ax.add_patch(Rectangle((0.07, VTA_INVERSION_CUTOFF), 6.2, 104 - VTA_INVERSION_CUTOFF,
                           color="#D6D6D6", alpha=0.5, zorder=0))
    ax.axhline(VTA_INVERSION_CUTOFF, color=INK, lw=2.4, ls="-.", zorder=1)
    for x, ls in ((1.0, "--"), (5.0, ":")):
        ax.axvline(x, color=INK, ls=ls, lw=2.4, zorder=1)
        ax.text(x, 108.5, f"{x:.0f}%", fontsize=15, fontweight="bold", color=INK,
                va="center", ha="center",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=INK, lw=1.2))

    for gene, x, y, reason, dy, ha, xm in T8:
        kept = reason is None
        ax.scatter([x], [y], s=540 if kept else 300,
                   marker="o", c=INK if kept else "white",
                   edgecolors=INK, linewidths=2.4, zorder=4)
        up = dy > 0
        ax.text(x * xm, y + dy, gene, ha=ha, va="bottom" if up else "top",
                fontsize=17, fontweight="bold", color=INK, zorder=5)
        if reason:
            ax.text(x * xm, y + dy + 4.0, reason, ha=ha, va="bottom",
                    fontsize=12.5, style="italic", color=MUTED, zorder=5,
                    linespacing=1.25)

    ax.set_xscale("log")
    ax.set_xlim(0.07, 6.2)
    ax.set_ylim(-5, 112)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.set_yticks([0, 25, 50, 75, 90, 100])
    ax.tick_params(labelsize=15)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontweight("bold")
    ax.set_xlabel("Worst SNc percentile  (lower = more central in SNc)",
                  fontsize=17, fontweight="bold", labelpad=12)
    ax.set_ylabel("Best VTA percentile  (higher = more peripheral in VTA)",
                  fontsize=17, fontweight="bold", labelpad=12)
    ax.grid(color="#EAEAEA", lw=1.0, zorder=0)
    ax.set_axisbelow(True)
    bare(ax)

    ax.text(0.078, 104.5, "KEPT: central in SNc, inverted in VTA",
            fontsize=17, fontweight="bold", color=INK, va="bottom")
    ax.set_title("Only four genes are both central in SNc and inverted in VTA",
                 fontsize=21, fontweight="bold", pad=18)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=INK, mec=INK, ms=18, label="Kept (4)"),
        Line2D([], [], marker="o", ls="", mfc="white", mec=INK, mew=2.4, ms=14,
               label="Rejected (9)"),
    ], loc="upper left", bbox_to_anchor=(0.012, 0.62), frameon=True,
        prop={"size": 16, "weight": "bold"})
    fig.tight_layout()
    save(fig, out_dir, "table8_quadrant")


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "C:/projects/rise_hdwgnca_data/table_figures"
    print(f"writing to {out_dir}")
    fig_table6(out_dir)
    fig_table7(out_dir)
    fig_table8(out_dir)


if __name__ == "__main__":
    main()
