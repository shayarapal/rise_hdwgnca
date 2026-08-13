# Mouse GSE233866 — CACNA1D Co-expression, SNc vs. VTA and Healthy vs. Lesioned

GSE233866 (Yaghmaeian Salmani et al., Karolinska Institute) profiles mouse midbrain dopaminergic
neurons across two arms: an **untreated** baseline cohort (6 animals, no disease) and a
**6-OHDA lesion** cohort (a different set of 6 animals, each contributing both a lesioned and an
intact/contralateral hemisphere). SNc/VTA labels are not provided by the dataset for either arm —
they're derived here via Sox6 (SNc) / Calb1 (VTA) marker scoring, mirroring the source paper's
own method (see `DATA_GSE233866/build_seurat_GSE233866.R` / `build_seurat_lesion_intact.R`).

Three sub-analyses, each answering a different question. See `analysis_results/SYNTHESIS.md`
for the full cross-study picture.

1. **`healthy_baseline/`** — SNc vs. VTA, no disease. *Are the regions structurally different at
   baseline?*
2. **`lesion_model/*_network` + `lesioned_vs_intact_combined_dme_snc/`** — lesioned vs. intact,
   within SNc. *Does disease change this region's network?*
3. **`lesion_model/vta_lesioned_network/` + `snc_vs_vta_combined_dme_lesioned/`** — SNc vs. VTA,
   within the lesioned state only. *What's specifically active in the vulnerable region during
   neurodegeneration itself* — the most direct analog of a region-selective treatment target.

---

## 1. `healthy_baseline/` — SNc vs. VTA, no disease

Population sizes: **2,226 SNc cells, 832 VTA cells** (pooled across all 6 untreated animals), out
of 4,514 total dopaminergic neurons identified via Th/Slc6a3/Ddc/Slc18a2 scoring. Soft power = 5
for both (independently determined per network, not assumed from the human dataset).

### Result: CACNA1D's module differs by region; connectivity is moderate in both, stronger in SNc

**Correction (2026-08-13):** this section previously reported kME=0.048 (SNc) and kME=0.094
(VTA), each "very weak." Both were misread from `modules.csv` — pulled from the wrong `kME_*`
column rather than each gene's actual assigned-module column. The correct values, re-verified
directly against `modules.csv` in both folders, are below.

| Region | Module | kME (own module) | Module size (genes) | # modules total |
|---|---|---|---|---|
| SNc | turquoise | 0.442 (moderate) | 1,352 | 9 |
| VTA | brown | 0.296 (weak-to-moderate) | 316 | 39 |

CACNA1D's module membership is a real, non-trivial connection in **both** mouse regions, not
noise — and notably, SNc's connectivity (0.44) is stronger than VTA's (0.30), the opposite
ranking of what this section previously (incorrectly) reported.

This baseline SNc value (0.44) also happens to be numerically higher than every kME in the
lesion-model section below (0.28–0.40) — **do not read that as a "connectivity drops once any
disease modeling happens" trend.** The healthy-baseline cohort (this section) and the lesion-arm
cohort (below) are two entirely separate sets of 6 animals from two different study arms; their
networks were built and clustered independently, from different cell counts, with different
module compositions. They are not a before/after series. The only *paired, controlled*
comparison in this project is intact-vs-lesioned **within** the lesion-arm cohort (0.346→0.398,
same animals, two hemispheres each) — that comparison, and the "connectivity increases after
lesioning" claim built on it, is correct and unaffected by this fix. The human PD-vs-control
comparison (kME 0.15 → 0.54) elsewhere in this project was independently re-checked against its
own `modules.csv` files as part of this fix and is correct as stated — this misread was isolated
to the two mouse healthy-baseline files. Treat the SNc-vs-VTA module identity shift
as informative on its own; don't over-read the exact kME magnitudes given the noise inherent to
per-gene estimates at this cell count. Similarly, the stark module-count difference (SNc: 9 vs.
VTA: 39, many under 100 genes) is likely a sample-size artifact (832 vs. 2,226 cells), not a
confirmed biological difference in network organization.

**Enrichment:** SNc's turquoise module hits "Dopaminergic synapse" (KEGG, p.adj=2.8×10⁻⁴)
directly; VTA's brown module hits "calcium ion transmembrane import into cytosol" (GO,
p.adj=7.8×10⁻³), on-topic for a calcium channel gene but a weaker signal.

### Formal preservation test (`snc_vs_vta_preservation/`, SNc as reference, 200 permutations)

Ran to completion easily (peaked at 3.88GB, ~12 minutes) — the memory ceiling that blocked this
same test on the human data was never a factor here. By WGCNA convention (Langfelder et al.
2011): Zsummary > 10 = strongly preserved, 2–10 = weak/moderate, < 2 = not preserved.

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

(`gold`/`grey` are WGCNA's internal random-control/unassigned-gene bins, excluded from ranking.
`moduleSize.obs` in `preservation.csv` is capped at 1000 genes — a standard WGCNA
computational-tractability limit, not each module's true size.)

**CACNA1D's module (turquoise) is strongly preserved (Zsummary=20.4), ranking 3rd of 8.** This
refines the qualitative finding above: even though CACNA1D gets a *different-colored* label when
VTA is analyzed independently (brown), the broader gene community it belongs to in SNc clearly
still holds together within VTA data — the module-identity shift likely reflects CACNA1D's own
relative connectivity shifting, not that neighborhood dissolving between regions.

All 6 animals cleared the `min_cells=100` metacell threshold, though VTA's smallest sample (s450,
102 cells) was close to that floor.

**Files** — `healthy_baseline/snc_network/`, `vta_network/`: `dendrogram.png`, `network_plot.png`
(kME plot), `module_eigengene_umap.png`, `module_scores_umap.png` (UCell),
`module_correlogram.png`/`.pdf`, `module_by_animal_heatmap.png` (mean hME per module per animal,
z-scored — spots animal outliers), `modules.csv`, `enrichr_table.csv` + `enrichr_plots/*.pdf`,
soft-power diagnostics. `healthy_baseline/snc_vs_vta_preservation/`: `preservation.csv`,
`preservation_plot.png`, `donor_counts.csv`, plus the reference (SNc) network's own outputs.

---

## 2. Lesioned vs. intact within SNc — the mouse PD-model comparison

Uses the study's disease arm: 6-OHDA-lesioned vs. contralateral intact hemisphere, same 6
animals (paired design) — the mouse analog of the human PD-vs-control axis.

Sanity check before any network analysis: 6-OHDA selectively kills dopaminergic neurons, and the
data confirms it — DA neuron counts dropped from 19,662 (intact) to 6,243 (lesioned) overall, and
**SNc lost proportionally more than VTA** (SNc: 8,637→2,291, 3.8x; VTA: 4,426→1,720, 2.6x),
matching this paper's whole premise of differential regional vulnerability.

### Two independently-built networks (`snc_intact_network/`, `snc_lesioned_network/`)

| Condition | n cells | Module | kME (own module) | Peak RAM |
|---|---|---|---|---|
| Intact | 8,637 | blue | 0.346 | 9.98GB |
| Lesioned | 2,291 | brown | 0.398 | 3.82GB |

CACNA1D's connectivity increases modestly under lesioning (0.346→0.398) — same *direction* as
the human PD result, though smaller. Enrichment is similar between the two (axon guidance,
synaptic transmission, cell adhesion) — less differentiated than the human data.

### Formal DME test (`lesioned_vs_intact_combined_dme_snc/`)

One combined network (both conditions pooled, 10,928 cells, 7.67GB peak), then `FindDMEs`.
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
statistically extreme than the equivalent human DME result (`human_GSE243639/
pd_vs_control_combined_dme/`, where CACNA1D's module does NOT reach significance). Combined with
the kME increase above: CACNA1D's broader co-expression program is suppressed in surviving
lesioned neurons, but CACNA1D itself becomes *relatively* more central to whatever remains of
that shrinking program — a hypothesis worth stating, not an established mechanism.

The brown module's log2FC is numerically unstable (module eigengenes can be negative/near-zero) —
trust the p-value, not that magnitude. Same caveat applies throughout every DME table in this
project.

---

## 3. SNc vs. VTA, within the lesioned state — what's active in the vulnerable region *during* disease

The most disease-and-region-specific comparison: not baseline regional structure (§1), not
disease-vs-healthy within one region (§2), but which genes/modules are more active in SNc than
VTA specifically while neurodegeneration is happening — the signal most directly relevant to a
region-selective treatment target.

**`lesion_model/vta_lesioned_network/`** (qualitative, mirrors §1/§2's separately-built
networks): CACNA1D lands in module blue, kME=0.283, 1,720 cells, 4.14GB peak.

**`lesion_model/snc_vs_vta_combined_dme_lesioned/`** (formal DME): one combined network pooling
SNc-lesioned (2,291 cells) + VTA-lesioned (1,720 cells) = 4,011 cells, 4.12GB peak. CACNA1D's
module here: **blue, kME=0.32**.

`FindDMEs(barcodes1=SNc, barcodes2=VTA)` — positive log2FC means higher in SNc:

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| greenyellow | −5.34 | ~0 | higher in VTA |
| turquoise | −6.70 | 7.2×10⁻²⁸⁸ | higher in VTA |
| tan | +4.70 | 1.1×10⁻²⁷¹ | higher in SNc |
| magenta | +6.83 | 3.4×10⁻¹⁴² | higher in SNc |
| **blue (CACNA1D)** | **+13.34** | **4.5×10⁻¹²³** | **higher in SNc, extremely significant** |
| pink | +5.41 | 2.0×10⁻¹¹⁹ | higher in SNc |
| green | +4.45 | 1.0×10⁻⁵¹ | higher in SNc |
| purple | −0.90 | 4.0×10⁻³² | higher in VTA |
| black | −1.36 | 3.4×10⁻³⁰ | higher in VTA |
| yellow | +12.78 | 8.8×10⁻²⁵ | higher in SNc |
| red | +0.28 | 1.5×10⁻¹⁷ | higher in VTA |
| brown | −9.43 | 0.068 (n.s.) | — |

**11 of 12 modules differ significantly between SNc and VTA during lesioning** — expected, since
these are fundamentally different neuron populations, not just disease-state variants of the
same one. What matters is that **CACNA1D's own module is among the strongest signals in the
entire table** (p.adj=4.5×10⁻¹²³) and is specifically *more active in SNc than VTA* during the
neurodegenerative process — the most direct answer this project has produced to "what's driving
SNc's selective vulnerability while it's actually happening." As with every DME table here,
treat the log2FC magnitudes as directional evidence, not literal fold-changes (module eigengenes
can be negative/near-zero, making that math numerically unstable) — the p-values are what to
trust.

---

## Caveats (apply throughout this folder)

- Sox6/Calb1-only cells were left unlabeled rather than forced into a region (1,456 of 4,514 DA
  neurons in the untreated arm) — a stricter or looser labeling rule would shift population sizes.
- Enrichr's `GO_Biological_Process_2023` / `KEGG_2021_Human` libraries are queried directly with
  mouse gene symbols (not converted to human orthologs); matching is generally case-insensitive
  in practice but hasn't been independently verified for these specific runs.
- These are mouse findings. Translating them back to the human results in
  `analysis_results/human_GSE243639/` is a separate, later inference — not something to claim
  directly from any comparison in this folder.
