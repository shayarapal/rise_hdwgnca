# DME: SNc vs. VTA, within Intact Cells — Statistics

Source: `dme_results.csv`, `soft_power_recommendation.json`. Sign convention (`DATA_GSE233866/run_dme_region_intact.R`, mirrors `run_dme_region_lesioned.R`): `barcodes1=SNc`, `barcodes2=VTA` → positive `avg_log2FC` = higher in SNc.

Built 2026-08-13 — fills the gap this project had: `snc_vs_vta_combined_dme_lesioned/` tested SNc vs VTA during lesioning, but nothing tested SNc vs VTA in the **intact** (unlesioned) hemisphere of the same lesion-arm cohort.

| Stat | Value |
|---|---|
| Pooled cells | 13,063 (8,637 SNc + 4,426 VTA) |
| Soft power (recommended) | 4 (`max_sft_r_sq=0.9669`) |
| Peak RAM | 10.17 GB |
| CACNA1D module | turquoise |
| CACNA1D kME (this network) | 0.327 |
| Module size | 1,629 genes |

## Differential module eigengenes

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| **turquoise (CACNA1D)** | **+13.21** | **≈0 (underflows double precision)** | higher in SNc |
| black | +7.69 | ≈0 | higher in SNc |
| yellow | −5.95 | ≈0 | higher in VTA |
| green | −3.24 | ≈0 | higher in VTA |
| red | +11.75 | 2.2×10⁻¹³³ | higher in SNc |
| brown | +1.40 | ≈0 | higher in SNc |
| pink | +0.50 | 6.9×10⁻³⁰ | higher in SNc |
| blue | +24.38 | 1.0 (n.s., unstable — see caveat) | — |

## Direct comparison: intact vs. lesioned SNc-vs-VTA

| | Intact (this file) | Lesioned (`../snc_vs_vta_combined_dme_lesioned/`) |
|---|---|---|
| CACNA1D module avg_log2FC | +13.21 | +13.34 |
| CACNA1D module p.adj | ≈0 (< double precision) | 4.5×10⁻¹²³ |
| Direction | higher in SNc | higher in SNc (same) |
| Pooled cells | 13,063 | 4,011 |

Essentially identical magnitude and direction in both states — SNc's version of this program
is not something that specifically emerges during lesioning; it is already present, at similar
strength, in the intact hemisphere.

## Gene-set overlap with related SNc networks

| Compared to | Shared genes | % of this module (1,629) |
|---|---|---|
| `../snc_vs_vta_combined_dme_lesioned/` (blue, its own module) | 823 | 51% |
| `../snc_intact_network/` (blue, 1,640 genes) | 1,198 | 74% |
| `../../healthy_baseline/snc_network/` (turquoise, 1,352 genes) | 853 | 52% |

Consistent with every other SNc network in this project — high, reproducible overlap across
independently-built networks from different cohorts and conditions.
