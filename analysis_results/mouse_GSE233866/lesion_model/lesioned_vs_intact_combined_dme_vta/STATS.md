# DME: Lesioned vs. Intact, within VTA — Statistics

Source: `dme_results.csv`, `soft_power_recommendation.json`. Sign convention (`DATA_GSE233866/run_dme_vta.R`, mirrors `run_dme_snc.R`): `barcodes1=lesioned`, `barcodes2=intact` → positive `avg_log2FC` = higher in lesioned.

Built 2026-08-13 — the VTA counterpart to `../lesioned_vs_intact_combined_dme_snc/`.

| Stat | Value |
|---|---|
| Pooled cells | 6,146 (4,426 intact + 1,720 lesioned) |
| Soft power (recommended) | 4 (`max_sft_r_sq` see `soft_power_recommendation.json`) |
| Peak RAM | 5.50 GB |
| CACNA1D module | green |
| CACNA1D kME (this network) | 0.265 |
| Module size | 191 genes |

## Differential module eigengenes (full table, 14 modules)

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| brown | −2.33 | 5.4×10⁻¹⁰⁸ | ↓ lesioned |
| purple | −1.70 | 1.8×10⁻⁴⁸ | ↓ lesioned |
| blue | +0.84 | 1.4×10⁻⁴⁰ | ↑ lesioned |
| magenta | −1.61 | 1.8×10⁻³⁹ | ↓ lesioned |
| pink | −5.68 | 2.0×10⁻²⁹ | ↓ lesioned |
| greenyellow | −0.60 | 2.2×10⁻²⁷ | ↓ lesioned |
| turquoise | +1.34 | 3.9×10⁻²² | ↑ lesioned |
| black | +0.95 | 5.3×10⁻¹⁵ | ↑ lesioned |
| tan | −0.65 | 9.3×10⁻⁷ | ↓ lesioned |
| **green (CACNA1D)** | **−1.23** | **2.2×10⁻⁵** | ↓ lesioned |
| yellow | +0.35 | 7.1×10⁻⁴ | ↑ lesioned |
| red | −0.17 | 5.7×10⁻³ | ↓ lesioned |
| cyan | +0.35 | 0.061 (n.s.) | — |
| salmon | −0.27 | 0.278 (n.s.) | — |

12 of 14 modules reach significance. CACNA1D's module (green) is significant but ranks
10th of 12 by p.adj — a real effect, not one of the table's strongest.

## Direct comparison to the SNc equivalent

| | SNc (`../lesioned_vs_intact_combined_dme_snc/`) | VTA (this file) |
|---|---|---|
| CACNA1D module avg_log2FC | −1.66 | −1.23 |
| CACNA1D module p.adj | 6.1×10⁻²⁶⁹ | 2.2×10⁻⁵ |
| Direction | ↓ lesioned | ↓ lesioned (same) |
| Rank among tested modules | 2nd most extreme p.adj of 7 | 10th of 12 |
| Pooled cells | 10,928 | 6,146 |

Same direction in both regions. SNc's effect is ~10²⁶⁴-fold more extreme by p.adj and
numerically larger by log2FC — not explainable by cell count alone (SNc has <2× VTA's
pooled cells here).

## Standalone-network kME comparison (descriptive, not a statistical test — see caveat in INTERPRETATION.tex)

| | SNc | VTA |
|---|---|---|
| Intact kME | 0.346 | 0.230 |
| Lesioned kME | 0.398 | 0.283 |
| Change | +0.052 (+15%) | +0.053 (+23%) |

Both regions increase by a similar absolute amount; VTA's relative increase is slightly larger.
