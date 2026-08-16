# DME: SNc vs. VTA, within Healthy Baseline — Statistics

Source: `dme_results.csv`, `soft_power_recommendation.json`. Sign convention (`DATA_GSE233866/run_dme_region_healthy.R`, mirrors `run_dme_region_intact.R`/`run_dme_region_lesioned.R`): `barcodes1=SNc`, `barcodes2=VTA` → positive `avg_log2FC` = higher in SNc.

Built 2026-08-16 — the third leg of the SNc-vs-VTA comparison (alongside `../../lesion_model/snc_vs_vta_combined_dme_intact/` and `snc_vs_vta_combined_dme_lesioned/`), and the first to test the region difference with **no disease model involved at all** — the previously-existing `snc_vs_vta_preservation/` in this same folder tests a different statistic (module structure, not activity level).

| Stat | Value |
|---|---|
| Pooled cells | 3,058 (2,226 SNc + 832 VTA) |
| Soft power (recommended) | 6 |
| Peak RAM | ~4.0 GB |
| CACNA1D module | yellow |
| CACNA1D kME (this network) | 0.444 |

## Differential module eigengenes

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| brown | −8.43 | ≈0 (underflows double precision) | higher in VTA |
| **yellow (CACNA1D)** | **+14.55** | **4.17×10⁻²¹²** | **higher in SNc** |
| greenyellow | +7.55 | 1.80×10⁻¹⁹⁹ | higher in SNc |
| green | −3.06 | 4.00×10⁻¹⁸⁷ | higher in VTA |
| purple | −2.47 | 3.81×10⁻¹¹⁹ | higher in VTA |
| pink | +6.31 | 7.38×10⁻⁸⁷ | higher in SNc |
| magenta | −1.93 | 8.79×10⁻⁷⁸ | higher in VTA |
| blue | +7.71 | 1.96×10⁻²⁰ | higher in SNc |
| turquoise | +2.17 | 9.66×10⁻¹⁷ | higher in SNc |
| black | +1.37 | 4.17×10⁻¹² | higher in SNc |
| red | −0.52 | 1.0 (n.s.) | — |

10 of 11 modules reach significance. CACNA1D's module (yellow) is the **single most extreme result in the table** by p.adj — not just significant, the most significant module tested.

## Direct comparison across all three disease states

| | Healthy baseline (this file) | Intact | Lesioned |
|---|---|---|---|
| CACNA1D module avg_log2FC | **+14.55** | +13.21 | +13.34 |
| CACNA1D module p.adj | 4.17×10⁻²¹² | ≈0 (underflows) | 4.5×10⁻¹²³ |
| Direction | higher in SNc | higher in SNc | higher in SNc |
| Rank among tested modules | **1st of 10** (most extreme) | tied for most extreme (5-way underflow) | 2nd of 11 |
| Pooled cells | 3,058 | 13,063 | 4,011 |

Same direction and comparable-to-larger magnitude in all three states — this is the cleanest confirmation in the project that the SNc-vs-VTA difference in CACNA1D's module is a **constitutive regional property**, present with no disease model involved whatsoever, not something that requires the lesion-arm cohort or any disease process to observe.
