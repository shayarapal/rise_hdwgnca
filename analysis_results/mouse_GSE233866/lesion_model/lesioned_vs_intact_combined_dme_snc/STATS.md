# DME: Lesioned vs. Intact, within SNc — Statistics

Source: `dme_results.csv`, `soft_power_recommendation.json`. Sign convention (`DATA_GSE233866/run_dme_snc.R`): `barcodes1=lesioned`, `barcodes2=intact` → positive `avg_log2FC` = higher in lesioned.

| Stat | Value |
|---|---|
| Pooled cells | 10,928 |
| Soft power (recommended) | 6 (`max_sft_r_sq=0.9676`) |
| Peak RAM | 7.67 GB |
| CACNA1D module | blue |
| CACNA1D kME (this network) | 0.40 |

## Differential module eigengenes

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| yellow | +6.93 | ~0 | ↑↑ lesioned |
| **blue (CACNA1D)** | **−1.66** | **6.1×10⁻²⁶⁹** | ↓ lesioned |
| red | −1.85 | 1.1×10⁻¹⁴¹ | ↓ lesioned |
| brown | unstable* | 7.1×10⁻¹⁷ | ↓ lesioned |
| turquoise | +1.82 | 2.1×10⁻¹⁶ | ↑ lesioned |
| black | −2.79 | 2.2×10⁻⁴ | ↓ lesioned |
| green | +0.75 | 1.0 (n.s.) | no change |

\* `avg_log2FC` numerically unstable near a zero-centered eigengene; `p.adj` remains trustworthy.
