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
`module_by_animal_heatmap.png` (mean harmonized module eigengene per module per animal, z-scored
per module — a quick way to spot animal-to-animal consistency or outliers within a region),
`modules.csv`, `enrichr_table.csv` + `enrichr_plots/*.pdf`, and soft-power diagnostics. All PNGs
are rendered large (2000-2800px) given VTA alone resolves to 39 modules.

## Formal preservation test (SNc as reference, VTA as query, 200 permutations)

Ran to completion easily on this dataset (peaked at 3.88GB, ~12 minutes) — the memory ceiling
that blocked this same test on the human data was never a factor here.

By WGCNA convention (Langfelder et al. 2011): Zsummary > 10 = strongly preserved, 2–10 =
weak/moderate, < 2 = not preserved.

| Module | Zsummary (preservation) | Rank (of 8 real modules) |
|---|---|---|
| brown | 53.1 | 1 |
| green | 47.2 | 2 |
| **turquoise (CACNA1D)** | **20.4** | **3** |
| yellow | 14.9 | 4 |
| blue | 14.6 | 5 |
| red | 10.9 | 6 |
| black | 7.9 | 7 |
| pink | 7.6 | 8 |

(`gold` and `grey` are WGCNA's internal random-control and unassigned-gene bins respectively —
not real modules, excluded from ranking. `moduleSize.obs` in `preservation.csv` is capped at
1000 genes for the largest modules — a standard WGCNA computational-tractability limit for the
permutation statistics, not each module's true gene count.)

**CACNA1D's module (turquoise) is strongly preserved (Zsummary=20.4), ranking 3rd of 8.** This
refines, rather than contradicts, the qualitative finding above: even though CACNA1D was
independently assigned to a *different-colored* module when VTA alone was analyzed (brown), the
broader gene community it belongs to in SNc clearly still holds together as a coherent
co-expression unit within VTA data. The module-identity shift likely reflects CACNA1D's own
*relative* connectivity/position shifting within an otherwise-stable neighborhood, rather than
that neighborhood itself dissolving between regions.

All 6 animals cleared the `min_cells=100` metacell threshold in both conditions, though VTA's
smallest sample (s450, 102 cells) was close to that floor — worth keeping in mind if this result
is revisited with more data.

Files: `snc_vs_vta_preservation/preservation.csv`, `preservation_plot.png`, plus the reference
(SNc) network's own `modules.csv`, `network_plot.png`, soft-power diagnostics, and
`donor_counts.csv`.

## Lesioned vs. Intact — the mouse PD-model comparison (disease-relevant)

Everything above uses the **untreated** (healthy baseline) arm of GSE233866 — no disease at all,
purely regional structure. This section uses the study's other arm: 6-OHDA-lesioned vs. the
contralateral intact hemisphere, a standard chemical PD model, in the same 6 animals (paired
design: each animal contributes both conditions). This is the actual mouse disease-vs-healthy
comparison, directly analogous to the human PD-vs-control axis.

Sanity check before any network analysis: 6-OHDA selectively kills dopaminergic neurons, and the
data confirms it — DA neuron counts dropped from 19,662 (intact) to 6,243 (lesioned) overall, and
**SNc lost proportionally more than VTA** (SNc: 8,637→2,291, 3.8x; VTA: 4,426→1,720, 2.6x),
matching this paper's whole premise of differential regional vulnerability.

### Two independently-built networks: SNc lesioned vs. SNc intact

| Condition | n cells | Module | kME (own module) | Peak RAM |
|---|---|---|---|---|
| Intact | 8,637 | blue | 0.346 | 9.98GB |
| Lesioned | 2,291 | brown | 0.398 | 3.82GB |

CACNA1D's connectivity increases modestly under lesioning (0.346→0.398) — same *direction* as
the human PD result, though a smaller shift. Enrichment is similar between the two (axon
guidance, synaptic transmission, cell adhesion) — less differentiated than the human data, which
showed a specific calcium-channel-signaling shift. Files: `snc_intact_network/`,
`snc_lesioned_network/` (same file set as the healthy `snc_network/`/`vta_network/` folders).

### Formal DME test: lesioned vs. intact, one combined SNc network

As with the human data, two independently-built networks aren't module-for-module comparable —
`snc_dme/` instead builds ONE combined network (both conditions pooled, 10,928 cells, 7.67GB
peak) and runs `FindDMEs` within it.

CACNA1D's module here: **blue, kME=0.40**.

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| yellow | +6.93 | ~0 | ↑↑ lesioned |
| **blue (CACNA1D)** | **−1.66** | **6.1×10⁻²⁶⁹** | **↓ lesioned, extremely significant** |
| red | −1.85 | 1.1×10⁻¹⁴¹ | ↓ lesioned |
| brown | (unstable, see caveat) | 7.1×10⁻¹⁷ | ↓ lesioned |
| turquoise | +1.82 | 2.1×10⁻¹⁶ | ↑ lesioned |
| black | −2.79 | 2.2×10⁻⁴ | ↓ lesioned |
| green | +0.75 | 1.0 (n.s.) | no change |

**CACNA1D's module is massively downregulated after lesioning** (p.adj ≈ 6×10⁻²⁶⁹) — far more
statistically extreme than the equivalent human PD-vs-control DME result (see
`analysis_results/human_snc_dme/`, where CACNA1D's module does NOT reach significance). Combined
with the kME increase above, this paints a coherent picture: CACNA1D's broader co-expression
program is suppressed in surviving lesioned neurons, but CACNA1D itself becomes *relatively* more
central to whatever remains of that shrinking program. Read as a hypothesis worth stating, not an
established mechanism — the top-level synthesis doc has the full cross-study interpretation.

The brown module's log2FC is numerically unstable (large magnitude, not a real fold-change) for
the same reason noted elsewhere: module eigengenes can be negative/near-zero. Trust the p-value,
not that number.

Files: `snc_dme/dendrogram.png`, `network_plot.png`, `modules.csv`, `dme_results.csv`,
`dme_volcano.png`, `dme_lollipop.png`, soft-power diagnostics.

## Caveats

- Sox6/Calb1-only cells were left unlabeled (1,456 of 4,514 DA neurons) rather than forced into
  a region — a stricter or looser labeling rule would shift these population sizes.
- Enrichr's `GO_Biological_Process_2023` / `KEGG_2021_Human` libraries are queried directly with
  mouse gene symbols (not converted to human orthologs first); matching is generally
  case-insensitive in practice but this hasn't been independently verified for this run.
- This is a mouse finding. Translating it back to the human SNc results in
  `analysis_results/{control_network,pd_network}/` is a separate, later inference — not
  something to claim directly from this comparison.
