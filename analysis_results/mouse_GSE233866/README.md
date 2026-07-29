# CACNA1D Co-expression Networks — Mouse SNc vs. VTA (GSE233866)

Same-species, same-study comparison (unlike the human control-vs-PD analysis one level up):
GSE233866 (Yaghmaeian Salmani et al., Karolinska Institute) profiles both SNc and VTA mouse
dopaminergic neurons together — 6 untreated DatCre;TRAP mice, ages P90/P530. SNc/VTA labels are
not provided by the dataset; they were derived here via Sox6 (SNc) / Calb1 (VTA) marker scoring,
mirroring the source paper's own method (see `DATA_GSE233866/build_seurat_GSE233866.R`).

Population sizes: **2,226 SNc cells, 832 VTA cells** (pooled across all 6 animals), out of 4,514
total dopaminergic neurons identified via Th/Slc6a3/Ddc/Slc18a2 scoring. Soft power = 5 for both
(independently determined per network, not assumed from the human dataset).

## Result: CACNA1D's module differs by region, but both connections are weak

| Region | Module | kME (own module) | Module size (genes) | # modules total |
|---|---|---|---|---|
| SNc | turquoise | 0.048 (very weak) | 1,352 | 9 |
| VTA | brown | 0.094 (weak) | 316 | 39 |

Unlike the human PD-vs-control result (kME 0.15 → 0.54), CACNA1D's module membership is weak in
**both** mouse regions. This likely reflects the much smaller cell counts here (thousands, not
tens of thousands) making per-gene kME estimates noisier — treat the SNc-vs-VTA module identity
shift itself as informative, but don't read strong quantitative weight into these particular kME
values.

The stark difference in module count (SNc: 9 modules vs. VTA: 39, many under 100 genes) is
likely driven by VTA's smaller sample size (832 vs. 2,226 cells) rather than a confirmed
biological difference in network organization — worth re-checking if a larger VTA dataset
becomes available.

## Enrichment

- **SNc turquoise** (CACNA1D's module): axon guidance, synaptic transmission, and — notably —
  **"Dopaminergic synapse"** (KEGG, p.adj=2.8×10⁻⁴), a direct hit for this cell type.
- **VTA brown** (CACNA1D's module): synaptic transmission and **calcium ion transmembrane
  import into cytosol** (GO, p.adj=7.8×10⁻³) — directly on-topic for a calcium channel gene,
  though a weaker signal than the SNc dopaminergic-synapse hit.

## Files

`snc_network/` and `vta_network/` each contain: `dendrogram.png`, `network_plot.png` (kME plot),
`module_eigengene_umap.png`, `module_scores_umap.png` (UCell), `module_correlogram.png`/`.pdf`,
`modules.csv`, `enrichr_table.csv` + `enrichr_plots/*.pdf`, and soft-power diagnostics.

## Caveats

- Sox6/Calb1-only cells were left unlabeled (1,456 of 4,514 DA neurons) rather than forced into
  a region — a stricter or looser labeling rule would shift these population sizes.
- Enrichr's `GO_Biological_Process_2023` / `KEGG_2021_Human` libraries are queried directly with
  mouse gene symbols (not converted to human orthologs first); matching is generally
  case-insensitive in practice but this hasn't been independently verified for this run.
- This is a mouse finding. Translating it back to the human SNc results in
  `analysis_results/{control_network,pd_network}/` is a separate, later inference — not
  something to claim directly from this comparison.
