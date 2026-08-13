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

CACNA1D: module **brown**, kME = 0.296 (316 genes in module) — weak-to-moderate, and weaker than
its SNc kME (0.442) — the reverse of what this document previously (incorrectly) reported.

**Correction (2026-08-13):** this previously read kME = 0.094 ("stronger than its SNc kME") —
misread from the wrong `kME_*` column in `modules.csv`. 0.296 is the actual `kME_brown` value
for Cacna1d, re-verified directly against the source file. The SNc/VTA connectivity ranking is
now the opposite of what was originally stated.

## Enrichment (significant, p.adj < 0.05)

Regenerated 2026-08-13 with the corrected `KEGG_2019_Mouse` library for KEGG terms (GO terms
were unaffected by that bug — see `DATA_GSE233866/build_vta_network.R`).

| Module | Term | DB | p.adj |
|---|---|---|---|
| brown | Calcium ion transmembrane import into cytosol | GO_Biological_Process_2023 | 7.8×10⁻³ |
