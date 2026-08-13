# SNc Network — Statistics

Source files in this folder: `soft_power_table.csv`, `soft_power_recommendation.json`, `modules.csv`, `enrichr_table.csv`

## Soft power selection

`sft_threshold=0.8`, `recommended_power=5`, `max_sft_r_sq=0.9644`

| Power | SFT.R.sq | mean.k. |
|---|---|---|
| 4 | 0.684 | 933.3 |
| **5** | **0.801** | **512.8** |
| 6 | 0.870 | 287.1 |

## Network summary

2,226 cells (pooled, 6 animals) → 9 modules.

CACNA1D: module **turquoise**, kME = 0.442 (1,352 genes in module) — moderate membership.

**Correction (2026-08-13):** this previously read kME = 0.048 ("weak, peripheral") — misread
from the wrong `kME_*` column in `modules.csv`. 0.442 is the actual `kME_turquoise` value for
Cacna1d, re-verified directly against the source file.

## Enrichment (significant, p.adj < 0.05)

Regenerated 2026-08-13 with the corrected `KEGG_2019_Mouse` library (previously queried against
`KEGG_2021_Human`, the wrong species — see `DATA_GSE233866/build_snc_network.R`).

| Module | Term | DB | p.adj |
|---|---|---|---|
| turquoise | Dopaminergic synapse | KEGG_2019_Mouse | 2.77×10⁻⁴ |
| turquoise | Axon guidance | KEGG_2019_Mouse | 2.77×10⁻⁴ |
