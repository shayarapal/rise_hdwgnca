# Mouse GSE233866 — SNc/VTA Healthy Baseline: Statistics Summary

All numbers below are read directly from the CSV/JSON files in this folder and
`../lesion_model/` — none are recalculated or approximated. Source file is noted
under each table. See `README.md` (this folder) for narrative interpretation, and
the published "hdWGCNA Field Guide" artifact for how to read each statistic.

## Population sizes

| Region | Cells (pooled, 6 animals) | Modules found | Soft power |
|---|---|---|---|
| SNc | 2,226 | 9 | 5 |
| VTA | 832 | 39 | 5 |

Source: `README.md`, `snc_network/soft_power_recommendation.json`, `vta_network/soft_power_recommendation.json`

## Soft power selection (SNc)

Source: `snc_network/soft_power_table.csv`, `snc_network/soft_power_recommendation.json`

| Power | SFT.R.sq | mean.k. |
|---|---|---|
| 1 | 0.277 | 6292.3 |
| 2 | 0.000 | 3268.1 |
| 3 | 0.303 | 1730.2 |
| 4 | 0.684 | 933.3 |
| **5 (recommended)** | **0.801** | **512.8** |
| 6 | 0.870 | 287.1 |
| 7 | 0.904 | 163.9 |
| 8 | 0.922 | 95.4 |
| 9 | 0.930 | 56.7 |
| 10 | 0.929 | 34.4 |

`sft_threshold` = 0.8, `max_sft_r_sq` = 0.9644 (peaks at higher powers but 5 is lowest to clear threshold).

## CACNA1D module membership, healthy baseline

| Region | Module | kME (own module) | Module size (genes) |
|---|---|---|---|
| SNc | turquoise | 0.048 | 1,352 |
| VTA | brown | 0.094 | 316 |

Source: `snc_network/modules.csv`, `vta_network/modules.csv`

## Preservation test (SNc = reference, VTA = query, 200 permutations)

Source: `snc_vs_vta_preservation/preservation.csv`, column `Zsummary.pres`
(`grey`/`gold` excluded — internal WGCNA control bins, not real modules)

| Module | Zsummary | Rank |
|---|---|---|
| brown | 53.1 | 1/8 |
| green | 47.2 | 2/8 |
| **turquoise (CACNA1D)** | **20.4** | **3/8** |
| yellow | 14.9 | 4/8 |
| blue | 14.6 | 5/8 |
| red | 10.9 | 6/8 |
| black | 7.9 | 7/8 |
| pink | 7.6 | 8/8 |

Interpretation convention (Langfelder et al. 2011): Zsummary > 10 = strongly preserved, 2–10 = weak/moderate, < 2 = not preserved.

## Donor cell counts (metacell QC)

Source: `snc_vs_vta_preservation/donor_counts.csv` — all 6 animals cleared `min_cells=100`

| Sample | SNc cells | VTA cells |
|---|---|---|
| s370 | 443 | 172 |
| s372 | 393 | 166 |
| s449 | 375 | 121 |
| s358 | 373 | 162 |
| s454 | 340 | 109 |
| s450 | 302 | 102 |

## Enrichment — significant hits (adjusted p < 0.05)

Source: `snc_network/enrichr_table.csv`, `vta_network/enrichr_table.csv`

| Region | Module | Term | DB | p.adj |
|---|---|---|---|---|
| SNc | turquoise | Dopaminergic synapse | KEGG_2021_Human | 2.8×10⁻⁴ |
| VTA | brown | Calcium ion transmembrane import into cytosol | GO_Biological_Process_2023 | 7.8×10⁻³ |

(Non-significant example row for contrast — SNc yellow, "Melanocyte Differentiation," p.adj = 0.273 — included in the Field Guide artifact, not a real hit.)

---

## Lesion model — DME statistics (`../lesion_model/`)

### Lesioned vs. intact, within SNc

Source: `../lesion_model/lesioned_vs_intact_combined_dme_snc/dme_results.csv`
Sign convention (`run_dme_snc.R`): `barcodes1=lesioned`, `barcodes2=intact` → positive = higher in lesioned.
CACNA1D's module in this pooled network: **blue**.

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| yellow | +6.93 | ~0 | ↑↑ lesioned |
| **blue (CACNA1D)** | **−1.66** | **6.1×10⁻²⁶⁹** | ↓ lesioned |
| red | −1.85 | 1.1×10⁻¹⁴¹ | ↓ lesioned |
| brown | unstable* | 7.1×10⁻¹⁷ | ↓ lesioned |
| turquoise | +1.82 | 2.1×10⁻¹⁶ | ↑ lesioned |
| black | −2.79 | 2.2×10⁻⁴ | ↓ lesioned |
| green | +0.75 | 1.0 (n.s.) | no change |

### Standalone networks (qualitative, not the formal DME test above)

Source: `../lesion_model/snc_intact_network/`, `../lesion_model/snc_lesioned_network/`

| Condition | n cells | Module | kME (own module) | Peak RAM |
|---|---|---|---|---|
| Intact | 8,637 | blue | 0.346 | 9.98 GB |
| Lesioned | 2,291 | brown | 0.398 | 3.82 GB |

### SNc vs. VTA, within lesioned cells only

Source: `../lesion_model/snc_vs_vta_combined_dme_lesioned/dme_results.csv`
Sign convention (`run_dme_region_lesioned.R`): `barcodes1=SNc`, `barcodes2=VTA` → positive = higher in SNc.
Pooled network: 2,291 SNc + 1,720 VTA = 4,011 cells, 4.12 GB peak. CACNA1D's module: **blue, kME=0.32**.

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| greenyellow | −5.34 | ~0 | higher in VTA |
| turquoise | −6.70 | 7.2×10⁻²⁸⁸ | higher in VTA |
| tan | +4.70 | 1.1×10⁻²⁷¹ | higher in SNc |
| magenta | +6.83 | 3.4×10⁻¹⁴² | higher in SNc |
| **blue (CACNA1D)** | **+13.34** | **4.5×10⁻¹²³** | higher in SNc |
| pink | +5.41 | 2.0×10⁻¹¹⁹ | higher in SNc |
| green | +4.45 | 1.0×10⁻⁵¹ | higher in SNc |
| purple | −0.90 | 4.0×10⁻³² | higher in VTA |
| black | −1.36 | 3.4×10⁻³⁰ | higher in VTA |
| yellow | +12.78 | 8.8×10⁻²⁵ | higher in SNc |
| red | +0.28 | 1.5×10⁻¹⁷ | higher in VTA |
| brown | unstable* | 0.068 (n.s.) | — |

\* `avg_log2FC` is numerically unstable when a module's eigengene averages near zero in either group (eigengenes can be negative — log2 of a near-zero-denominator ratio is not a reliable magnitude). The `p.adj` value is unaffected and remains trustworthy.

### Lesion-driven cell loss (sanity check, not a module statistic)

Source: `README.md` (this folder)

| Region | Intact | Lesioned | Fold loss |
|---|---|---|---|
| SNc | 8,637 | 2,291 | 3.8× |
| VTA | 4,426 | 1,720 | 2.6× |
