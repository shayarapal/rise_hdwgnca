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

## Module preservation

| Field | Value |
|---|---|
| hdWGCNA function | `ModulePreservation(seurat_obj, seurat_ref, wgcna_name)` |
| New endpoints | `POST /module-preservation` → job_id; `GET /results/{job_id}/preservation` → JSON |
| Extra input | Path to a second reference `.h5Seurat` file |
| Notes | Computationally intensive — may need `workers = 4` in `entrypoint.R` |

---

## Consensus network (multi-group)

| Field | Value |
|---|---|
| hdWGCNA function | `ConstructNetwork(seurat_obj, consensus = TRUE, ...)` |
| Change to existing | Add `consensus` boolean to `/analyze` request body (default `false`) |
| Notes | No new endpoint needed — parameterise the existing one |
