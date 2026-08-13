# VTA Network — Statistics

Source files in this folder: `soft_power_table.csv`, `soft_power_recommendation.json`, `modules.csv`, `enrichr_table.csv`

## Soft power selection

`sft_threshold=0.8`, `recommended_power=5`, `max_sft_r_sq=0.9766`

| Power | SFT.R.sq | mean.k. |
|---|---|---|
| 4 | 0.719 | 874.3 |
| **5** | **0.850** | **482.6** |
| 6 | 0.895 | 271.8 |

## Network summary

832 cells (pooled, 6 animals) → 39 modules.

CACNA1D: module **brown**, kME = 0.094 (316 genes in module) — weak, but stronger than its SNc kME.

## Enrichment (significant, p.adj < 0.05)

| Module | Term | DB | p.adj |
|---|---|---|---|
| brown | Calcium ion transmembrane import into cytosol | GO_Biological_Process_2023 | 7.8×10⁻³ |
