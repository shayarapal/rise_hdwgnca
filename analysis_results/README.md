# CACNA1D Co-expression Networks — Control vs. PD (SNc Dopaminergic Neurons)

Data: GSE243639 (29 donors: 15 Parkinson's, 14 control; substantia nigra pars compacta snRNA-seq).
Pipeline: `r-service` hdWGCNA bridge (`SetupForWGCNA` → `MetacellsByGroups` → `NormalizeMetacells` →
`SetDatExpr` → `TestSoftPowers` → `ConstructNetwork` → `ModuleEigengenes` → `ModuleConnectivity`),
restricted to the "Dopaminergic Neurons" cell type, run separately per condition. Soft power = 5
for both (network_type="signed", gene_select="fraction" 0.05, k=25, max_shared=15).

## Result: CACNA1D's module assignment differs by condition

| Condition | Module   | kME (own module) | Module size (genes) |
|-----------|----------|-------------------|----------------------|
| Control   | turquoise | 0.147 (weak)      | 3,885                |
| PD        | blue      | 0.540 (strong)    | 1,511                |

CACNA1D is only weakly tied to its assigned module in control tissue, but is a strongly-connected
member of the "blue" module in PD tissue — a qualitative shift in its co-expression neighborhood
between conditions.

## Files

- `control_network/modules.csv` — full gene→module table, control condition
- `control_network/dendrogram.png` — gene clustering dendrogram with module color bar (`PlotDendrogram`)
- `control_network/network_plot_kME.png` — per-module kME plot (`PlotKMEs`), control condition
- `control_network/module_eigengene_umap.png` — harmonized module eigengene (hME) activity projected
  onto the single-cell UMAP, one panel per module (`ModuleFeaturePlot`); high values = cells where
  that module's genes are collectively highly active
- `control_network/module_correlogram.png` (+ `.pdf`) — pairwise module eigengene correlation
  heatmap (`ModuleCorrelogram`); modules that cluster together likely reflect related biology
- `control_network/module_scores_umap.png` — per-cell UCell module scores (top 25 hub genes per
  module, rank-based, outlier-robust) projected onto the UMAP (`ModuleExprScore` + `ModuleFeaturePlot`)
- `control_network/soft_power_plot.png`, `soft_power_table.csv`, `soft_power_recommendation.json`
- `pd_network/` — same eight outputs, PD condition

**hME vs. module score:** use hMEs (`module_eigengene_umap.png`) for statistical analyses —
they're continuous and directly comparable across cells. Use UCell scores
(`module_scores_umap.png`) for visualization — more intuitive, less sensitive to outlier cells.

## What this is and isn't

This is **two independently-built networks**, not a formal preservation test. It shows each
condition's own module structure and where CACNA1D falls in each — it does *not* give a
statistical answer to "is CACNA1D's module preserved between conditions" (that requires
`ModulePreservation()`'s permutation test, comparing the network built in one condition against
data from the other).

That formal test (`POST /module-preservation`, `ref_group=control`, `query_group=PD`) was attempted
repeatedly but did not complete on this machine: it requires holding both conditions' processed
networks in memory simultaneously through a permutation loop, and consistently exceeded the
available RAM (raised as high as 44GB) somewhere between permutation ~20 and ~80 out of 250,
non-deterministically. Two real memory bugs were found and fixed along the way (a hdWGCNA
sparse-matrix inefficiency in `SelectNetworkGenes`, patched in `r-service/seuratdisk_compat.R`;
and a `plumber.R` inefficiency holding both conditions' raw data in memory unnecessarily early) —
what remains is a genuine memory requirement of the permutation test itself on this dataset size,
not a bug in this project's code.
