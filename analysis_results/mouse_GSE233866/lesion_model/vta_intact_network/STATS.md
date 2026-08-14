# VTA, Intact Hemisphere — Statistics

Source: `soft_power_recommendation.json`, `modules.csv`, `enrichr_table.csv`. Built
2026-08-13 — the VTA counterpart to `../snc_intact_network/`, which previously had no
equivalent (see `../vta_lesioned_network/STATS.md` for context on why this gap existed).

| Stat | Value |
|---|---|
| Cells | 4,426 |
| Soft power (recommended) | 4 (`max_sft_r_sq=0.9593`) |
| CACNA1D module | pink |
| CACNA1D kME (own module) | 0.230 |
| Module size | 132 genes |
| Peak RAM | 4.23 GB |

## Enrichment

**No terms reach significance** (best: "Morphine addiction", KEGG_2019_Mouse, p.adj = 0.120;
826 terms tested total). This is the only CACNA1D-containing module in the entire project
with zero significant enrichment hits — every other one (SNc healthy/intact/lesioned, VTA
healthy/lesioned) has at least one term at p.adj < 0.02.

## Gene-set overlap with other VTA networks

| Compared to | Shared genes | % of this module (132) |
|---|---|---|
| `../vta_lesioned_network/` (blue, 695 genes) | 19 | 14% |
| `../lesioned_vs_intact_combined_dme_vta/` (green, 191 genes) | 84 | 64% |
| `../../healthy_baseline/vta_network/` (brown, 316 genes) | 40 | 30% |

For comparison, the equivalent SNc pair (`snc_intact_network` vs. `snc_lesioned_network`)
shares 792 of 1,640 genes (48%) — VTA's intact-vs-lesioned overlap (14%) is roughly a third
of that.
