#!/usr/bin/env python3
"""Regenerate the paper's DME volcano figures, eNeuro-compliant.

Reads the committed dme_results.csv tables and redraws them. Three eNeuro rules drive the
styling and are the reason this script exists rather than the ad-hoc originals:

  - "Remove top and right borderlines ... do not box the panels in"  -> left/bottom spines only
  - "Figures with red and green ... convert to magenta and green"    -> module colours dropped
  - two-bar graphs are not permitted                                 -> the old Figure 2 is gone

WGCNA module colours are arbitrary labels with no quantitative meaning, so encoding them as
point fills bought nothing and was what put red and green in the same panel. Points are now
neutral grey, CACNA1D's module is magenta, and the module name stays as the text label.

Shape still encodes plotting status, which does carry meaning:
  circle        p.adj plotted directly
  triangle-up   p.adj underflowed double precision (reported as 0), drawn at the ceiling
  triangle-left avg_log2FC beyond the axis, drawn at the edge

Usage:  python make_figures.py [RESULTS_DIR] [OUT_DIR]
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "pdf.fonttype": 42,          # embed TrueType; outline in a vector editor if EPS is required
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "axes.linewidth": 1.0,
})

GREY, GREY_EDGE = "#B0B0B0", "#4D4D4D"
MAGENTA, MAGENTA_EDGE = "#D81B9A", "#000000"

DEFAULT_RESULTS = os.environ.get(
    "ANALYSIS_RESULTS",
    "C:/projects/rise_hdwgnca_data/analysis_results/mouse_GSE233866",
)

# panel -> (results subdir, CACNA1D's module in that network, title, subtitle)
FIG1 = [
    ("healthy_baseline/snc_vs_vta_combined_dme_healthy", "yellow", "Healthy baseline",
     "SNc vs. VTA - untreated cohort, 3,058 pooled cells"),
    ("lesion_model/snc_vs_vta_combined_dme_intact", "turquoise", "Intact",
     "SNc vs. VTA - lesion-arm cohort, unlesioned side, 13,063 pooled cells"),
    ("lesion_model/snc_vs_vta_combined_dme_lesioned", "blue", "Lesioned",
     "SNc vs. VTA - lesion-arm cohort, 6-OHDA side, 4,011 pooled cells"),
]
FIG2 = [
    ("lesion_model/lesioned_vs_intact_combined_dme_snc", "blue", "SNc",
     "Lesioned vs. Intact, 10,928 pooled cells"),
    ("lesion_model/lesioned_vs_intact_combined_dme_vta", "green", "VTA",
     "Lesioned vs. Intact, 6,146 pooled cells"),
]


def read_dme(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [{"module": r["module"],
                 "lfc": float(r["avg_log2FC"]),
                 "padj": float(r["p_val_adj"])}
                for r in csv.DictReader(fh)]


def fmt_p(p):
    """Render p.adj the way the captions do: mantissa x 10^exp."""
    if p == 0:
        return "< 1x10$^{-300}$"
    e = math.floor(math.log10(p))
    return f"{p / 10 ** e:.2f}x10$^{{{e}}}$"


def repel(ax, items, xlim, ylim, obstacles=(), rounds=520):
    """Nudge labels apart in axes-fraction space, then draw a leader line to each point.

    `obstacles` are fixed (x0, y0, x1, y1) axes-fraction rectangles - the CACNA1D callout
    box - that labels are pushed clear of.
    """
    xr = xlim[1] - xlim[0]
    yr = ylim[1] - ylim[0]
    pts = [((it["x"] - xlim[0]) / xr, (it["y"] - ylim[0]) / yr) for it in items]
    # points near the ceiling get their label below, otherwise it lands in the title band
    labs = [[x, y - 0.058 if y > 0.86 else y + 0.058] for x, y in pts]
    w, h = 0.092, 0.045            # approximate label box, axes fractions

    for _ in range(rounds):
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                dx = labs[j][0] - labs[i][0]
                dy = labs[j][1] - labs[i][1]
                ox = w - abs(dx)
                oy = h - abs(dy)
                if ox > 0 and oy > 0:          # overlap -> push along the shallower axis
                    if ox < oy:
                        s = math.copysign(ox / 2, dx or 1.0)
                        labs[i][0] -= s
                        labs[j][0] += s
                    else:
                        s = math.copysign(oy / 2, dy or 1.0)
                        labs[i][1] -= s
                        labs[j][1] += s
            for x0, y0, x1, y1 in obstacles:   # evict from the callout box
                if x0 - w / 2 < labs[i][0] < x1 + w / 2 and y0 - h / 2 < labs[i][1] < y1 + h / 2:
                    left, right = labs[i][0] - (x0 - w / 2), (x1 + w / 2) - labs[i][0]
                    down, up = labs[i][1] - (y0 - h / 2), (y1 + h / 2) - labs[i][1]
                    m = min(left, right, down, up)
                    if m == left:
                        labs[i][0] = x0 - w / 2
                    elif m == right:
                        labs[i][0] = x1 + w / 2
                    elif m == down:
                        labs[i][1] = y0 - h / 2
                    else:
                        labs[i][1] = y1 + h / 2
            for k, (lo, hi) in enumerate(((0.012, 0.988), (0.012, 0.975))):
                labs[i][k] = min(max(labs[i][k], lo), hi)

    # Uncross leader lines: if trading two labels shortens the total run, they were crossed.
    def d2(a, b):
        return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2

    for _ in range(6):
        swapped = False
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                if d2(pts[i], labs[j]) + d2(pts[j], labs[i]) < d2(pts[i], labs[i]) + d2(pts[j], labs[j]):
                    labs[i], labs[j] = labs[j], labs[i]
                    swapped = True
        if not swapped:
            break

    # Final settle: keep labels off the markers themselves, and off each other again.
    for _ in range(140):
        for i, lab in enumerate(labs):
            for px, py in pts:
                dx, dy = lab[0] - px, lab[1] - py
                ox, oy = 0.048 - abs(dx), 0.034 - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy:
                        lab[0] += math.copysign(ox, dx or 1.0)
                    else:
                        lab[1] += math.copysign(oy, dy or 1.0)
            for j, other in enumerate(labs):
                if i == j:
                    continue
                dx, dy = lab[0] - other[0], lab[1] - other[1]
                ox, oy = w - abs(dx), h - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy:
                        lab[0] += math.copysign(ox / 2, dx or 1.0)
                    else:
                        lab[1] += math.copysign(oy / 2, dy or 1.0)
            lab[0] = min(max(lab[0], 0.012), 0.988)
            lab[1] = min(max(lab[1], 0.012), 0.975)

    for it, (px, py), (lx, ly) in zip(items, pts, labs):
        ax.plot([lx * xr + xlim[0], px * xr + xlim[0]],
                [ly * yr + ylim[0], py * yr + ylim[0]],
                color="#808080", lw=0.6, zorder=1)
        ax.text(lx * xr + xlim[0], ly * yr + ylim[0], it["module"],
                ha="center", va="center", zorder=5,
                fontsize=10.5 if it["hit"] else 9.5,
                fontweight="bold" if it["hit"] else "normal",
                color=MAGENTA_EDGE if it["hit"] else "#1A1A1A")


def panel(ax, rows, target, title, subtitle, xlim, ceiling, ylim, left_label, right_label):
    items = []
    for r in rows:
        under = r["padj"] == 0
        y = ceiling if under else -math.log10(r["padj"])
        off = r["lfc"] < xlim[0] or r["lfc"] > xlim[1]
        x = min(max(r["lfc"], xlim[0] + 0.12), xlim[1] - 0.12)
        hit = r["module"] == target
        ax.scatter([x], [y], marker="<" if off else ("^" if under else "o"),
                   s=235 if hit else 125,
                   c=MAGENTA if hit else GREY,
                   edgecolors=MAGENTA_EDGE if hit else GREY_EDGE,
                   linewidths=2.0 if hit else 0.9, zorder=4 if hit else 3)
        items.append({"module": r["module"], "x": x, "y": y, "hit": hit})
        if off:
            ax.annotate(f"log2FC={r['lfc']:.1f}\n(off-scale)", (x, y), xytext=(6, -17),
                        textcoords="offset points", fontsize=7.5, color="#707070")

    ax.axhline(-math.log10(0.05), ls="--", lw=0.9, color="#9A9A9A", zorder=0)
    ax.axvline(0, lw=0.9, color="#D0D0D0", zorder=0)

    # Lower right: the only region empty in all five panels, and clear of the ceiling row.
    tgt = next(r for r in rows if r["module"] == target)
    ax.text(0.985, 0.30,
            f"log$_2$FC = {tgt['lfc']:+.2f}\np.adj = {fmt_p(tgt['padj'])}\n{target} (CACNA1D)",
            transform=ax.transAxes, va="top", ha="right", fontsize=10.5,
            fontweight="bold", color=MAGENTA_EDGE,
            bbox=dict(boxstyle="round,pad=0.42", fc="white", ec=MAGENTA, lw=1.1, alpha=0.95))
    callout = (0.70, 0.13, 0.99, 0.31)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.spines["top"].set_visible(False)      # eNeuro: no top/right borderlines
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.set_title(title, fontsize=16, fontweight="bold", pad=22)
    ax.text(0.5, 1.012, subtitle, transform=ax.transAxes, ha="center",
            fontsize=9, color="#555555")
    ax.set_xlabel("Average log$_2$ (Fold Change)\n"
                  f"$\\leftarrow$ higher in {left_label}      "
                  f"higher in {right_label} $\\rightarrow$", fontsize=10.5)
    repel(ax, items, xlim, ylim, obstacles=(callout,))


def build(spec, results_dir, out_dir, name, suptitle, xlim, left_label, right_label):
    data = [(read_dme(os.path.join(results_dir, d, "dme_results.csv")), t, ti, s)
            for d, t, ti, s in spec]
    max_finite = max(-math.log10(r["padj"]) for rows, *_ in data for r in rows if r["padj"] > 0)
    ceiling = math.ceil(max_finite * 1.06 / 10) * 10
    ylim = (-ceiling * 0.055, ceiling * 1.16)

    fig, axes = plt.subplots(1, len(data), figsize=(7.6 * len(data), 8.4), squeeze=False)
    for ax, (rows, target, title, subtitle) in zip(axes[0], data):
        panel(ax, rows, target, title, subtitle, xlim, ceiling, ylim, left_label, right_label)
    axes[0][0].set_ylabel("$-$log$_{10}$ (Adjusted P-value)", fontsize=12)

    fig.suptitle(suptitle, fontsize=17, y=0.985)
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=GREY, mec=GREY_EDGE, ms=10,
               label="Module (p.adj plotted directly)"),
        Line2D([], [], marker="^", ls="", mfc=GREY, mec=GREY_EDGE, ms=10,
               label=f"p.adj underflowed, drawn at ceiling ({ceiling})"),
        Line2D([], [], marker="<", ls="", mfc=GREY, mec=GREY_EDGE, ms=10,
               label="log$_2$FC off-scale, drawn at edge"),
        Line2D([], [], marker="o", ls="", mfc=MAGENTA, mec=MAGENTA_EDGE, mew=2, ms=12,
               label="CACNA1D's module"),
    ], loc="upper center", bbox_to_anchor=(0.5, 0.945), ncol=4, frameon=True, fontsize=10)
    fig.text(0.5, 0.022, "p.adj = 0.05 (dashed line)", ha="center", fontsize=8.5, color="#8A8A8A")
    fig.tight_layout(rect=[0, 0.035, 1, 0.915])

    os.makedirs(out_dir, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), format=ext,
                    bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print(f"{name}: {len(data)} panels, ceiling={ceiling} -> {out_dir}/{name}.{{pdf,png}}")


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_RESULTS
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(results_dir, "figures")
    build(FIG1, results_dir, out_dir, "figure1",
          "CACNA1D's module is more active in SNc than VTA in every cohort, "
          "with and without a disease model",
          (-12.0, 27.0), "VTA", "SNc")
    build(FIG2, results_dir, out_dir, "figure2",
          "CACNA1D's module activity falls after lesioning in both regions - "
          "far more strongly in SNc",
          (-8.0, 8.0), "intact", "lesioned")


if __name__ == "__main__":
    main()
