# DME: SNc vs. VTA, within Lesioned Cells — Statistics

Source: `dme_results.csv`, `soft_power_recommendation.json`. Sign convention (`DATA_GSE233866/run_dme_region_lesioned.R`): `barcodes1=SNc`, `barcodes2=VTA` → positive `avg_log2FC` = higher in SNc.

| Stat | Value |
|---|---|
| Pooled cells | 4,011 (2,291 SNc + 1,720 VTA) |
| Soft power (recommended) | 5 (`max_sft_r_sq=0.9424`) |
| Peak RAM | 4.12 GB |
| CACNA1D module | blue |
| CACNA1D kME (this network) | 0.32 |

## Differential module eigengenes

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

\* `avg_log2FC` numerically unstable near a zero-centered eigengene; `p.adj` remains trustworthy (and here, not significant regardless).

11 of 12 modules differ significantly between SNc and VTA during lesioning.
