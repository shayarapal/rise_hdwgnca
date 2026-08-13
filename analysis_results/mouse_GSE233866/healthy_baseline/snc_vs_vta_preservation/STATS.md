# SNc → VTA Preservation — Statistics

Reference network: SNc. Query: VTA. 200 permutations. Source: `preservation.csv` (`Zsummary.pres` column), `donor_counts.csv`.

## Zsummary ranking (real modules only; `grey`/`gold` excluded as internal control bins)

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

Convention (Langfelder et al. 2011): Zsummary > 10 strongly preserved, 2–10 weak/moderate, < 2 not preserved.

## Donor QC (min_cells = 100 threshold)

| Sample | SNc cells | VTA cells |
|---|---|---|
| s370 | 443 | 172 |
| s372 | 393 | 166 |
| s449 | 375 | 121 |
| s358 | 373 | 162 |
| s454 | 340 | 109 |
| s450 | 302 | 102 |

Peak run: 3.88 GB, ~12 minutes.
