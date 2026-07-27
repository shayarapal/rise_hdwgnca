# Planned Skills (Future Endpoints)

This document catalogs hdWGCNA capabilities not yet wrapped as endpoints. Each entry
follows the same addition pattern as the existing endpoints: file-path in, `job_id`
out, poll `/status/{job_id}`, fetch from `/results/{job_id}/<result>`.

When implementing, use `PROMPT_TEMPLATE.md` and add a CHANGELOG entry.

---

## Hub gene export

| Field | Value |
|---|---|
| hdWGCNA function | `GetHubGenes(seurat_obj, n_hubs, wgcna_name)` |
| New endpoint | `GET /results/{job_id}/hub-genes` |
| Output | JSON array: `[{ gene, module_color, kME }]` |
| Notes | `n_hubs` (top N per module) should be a query parameter, default 10 |

---

## Trait correlation

| Field | Value |
|---|---|
| hdWGCNA functions | `ModuleTraitCorrelation()`, `ModuleTraitCorHeatmap()` |
| New endpoints | `POST /trait-correlation` → job_id; `GET /results/{job_id}/trait-heatmap` → PNG |
| Extra input | Path to phenotype/trait CSV (numeric columns = traits, rows = samples) |
| Notes | Requires `group.by.vars` to match the grouping used in `ModuleEigengenes` |

---

## Gene set enrichment

| Field | Value |
|---|---|
| hdWGCNA function | `OverlapModulesWith(seurat_obj, enrichment_db, wgcna_name)` |
| New endpoint | `GET /results/{job_id}/enrichment` → JSON |
| Extra input | Path to gene set `.gmt` file on shared volume |
| Notes | Confirm exact `OverlapModulesWith` signature before implementing |

---

## Module UMAP

| Field | Value |
|---|---|
| hdWGCNA functions | `RunModuleUMAP()`, `ModuleUMAPPlot()` |
| New endpoints | `POST /module-umap` → job_id; `GET /results/{job_id}/umap-plot` → PNG |
| Notes | Requires a completed `/analyze` job (Seurat object with modules). Consider persisting the Seurat RDS in `out_dir` from `/analyze` so it can be reloaded here. |

---

## Module preservation via NetRep (faster alternative)

| Field | Value |
|---|---|
| hdWGCNA function | `ModulePreservationNetRep(seurat_query, seurat_ref, name, n_permutations, n_threads, TOM_use, wgcna_name, wgcna_name_ref)` |
| Change to existing | A `method` switch on the existing `POST /module-preservation` |
| Extra input | None, but the **query** needs its own network (`ConstructNetwork`) so a query TOM exists for `TOM_use` |
| Notes | Needs the `NetRep` package, not currently in `install.R`. Plots via `PlotModulePreservationLollipop`. The permutation-based `ModulePreservation` already wrapped is the tutorial's default path; add this only if runtime becomes the bottleneck. |

---

## Consensus network (multi-group)

| Field | Value |
|---|---|
| hdWGCNA function | `ConstructNetwork(seurat_obj, consensus = TRUE, ...)` |
| Change to existing | Add `consensus` boolean to `/analyze` request body (default `false`) |
| Notes | No new endpoint needed — parameterise the existing one |
