# Changelog

All notable changes to this project will be documented in this file.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).  
Versioning follows [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Known issues
- `MetacellsByGroups` still uses `reduction = "pca"` and `max_shared = 15`, where the hdWGCNA
  tutorial uses `reduction = 'harmony'` and `max_shared = 10`. **Deliberate, not an oversight.**
  Metacells never span donors (`group_by` includes `sample_id`), so the uncorrected reduction has
  limited effect. Revisit if metacells are ever built across batches.
- `build_seurat_GSE243639.R` labels a mixed neuronal cluster as `Dopaminergic Neurons`
  (it wins the TH/SLC6A3 panel with a weak score of ~0.33 while only ~40 % of its nuclei
  express either marker). `relabel_cluster7.R` corrects this after the fact and **must be
  re-run after any rebuild**. A purity guard in the labelling step would fix it properly.
- Only **11 of 29 donors** (4 PD, 7 control) clear `MetacellsByGroups(min_cells = 100)` for
  dopaminergic neurons — PD brains have lost the very cells being counted. Fine for a pooled
  network, but a PD-vs-control comparison built on 4 PD donors is underpowered. A property of
  the data, not a bug.

---

## [0.4.0] - 2026-07-13

> **Network results change.** The co-expression network now spans 7,465 genes instead of 2,000.
> Any modules produced before this release were built on a different gene universe and should be
> regenerated.

### Fixed
- **`SetupForWGCNA` was not selecting genes the way the hdWGCNA tutorial does.** `plumber.R`
  omitted `gene_select` and `fraction`, so `SelectNetworkGenes` fell back to its
  `gene_select = "variable"` default — the 2,000 `FindVariableFeatures` genes, computed across
  *all 107,113 cells of every cell type*, i.e. genes that separate oligodendrocytes from
  astrocytes. The tutorial passes `gene_select = "fraction", fraction = 0.05` (genes expressed in
  ≥ 5 % of cells in the analysed group), which for the dopaminergic neurons is **7,465 genes**.
  hdWGCNA can only find a gene's partners among genes in the network, so the old setting asked
  "what co-expresses with CACNA1D?" while excluding ~5,500 genes dopaminergic neurons express.
  CACNA1D itself (8.0 % of DA nuclei) was present under both settings, so no prior result was
  *invalid* — just narrow. **6,826 genes are newly searchable as CACNA1D partners.**
- **`/analyze` built the network on different metacells than the soft power was chosen from.**
  `MetacellsByGroups` is kNN-based and takes no `seed` argument, and none was set; `/analyze`
  reloads the source object and recomputes from scratch. The soft power picked from the
  `/test-soft-powers` curve was therefore applied to a *different* metacell set, silently.
  `set.seed(42)` now precedes both `MetacellsByGroups` calls — verified to reproduce identical
  metacell barcodes (5,482 metacells) across runs.
- **TOM files escaped the job's output directory.** `ConstructNetwork` defaults to
  `tom_outdir = "TOM"` — a *relative* path — so the topological overlap matrix was written to the
  r-service working directory, not `out_dir`. With `tom_name = NULL`, `overwrite_tom = FALSE`, and
  `plan(multisession, workers = 2)`, a rerun could trip on a stale TOM and concurrent jobs could
  collide. Now passes `tom_outdir = p$out_dir` and `tom_name = p$wgcna_name`.
- `frontend/src/App.jsx` defaulted `group_by` to `['cell_type', 'Sample']`; `Sample` is not a
  column in this dataset. Now `['cell_type', 'sample_id']`, with the `ParamForm` placeholder to
  match.

### Removed
- `plumber.R` wrote `seurat_after_soft_powers.rds` with the comment *"Persist object so /analyze
  can reload it if needed"* — `/analyze` never read it. Dead code; with the metacell seed in place
  it is not needed.

---

## [0.3.0] - 2026-07-13

### Added
- `DATA_UNZIPPED/make_da_object.R` — writes `seurat_GSE243639_SNc_DA.rds`, a dopaminergic-only
  analysis object. Subsets to the 8,216 DA nuclei and drops the `scale.data` layer, taking the
  object from **7.75 GB to 0.33 GB in RAM** (23×) and 0.06 GB on disk. Retains the `pca`
  reduction (`MetacellsByGroups` hardcodes `reduction = "pca"`) and the 2,000 variable features
  (`SelectNetworkGenes` defaults to `gene_select = "variable"`). No hdWGCNA result changes.
- `.claude/skills/release/SKILL.md` — `/release` skill: bumps the version, writes the CHANGELOG
  entry, commits, tags, and pushes, with a hard guard against committing the ~11 GB dataset.

### Fixed
- **`vector memory limit of 32.0 Gb` when running the pipeline.** The r-service runs
  `plan(multisession, workers = 2)`, so two worker processes each held a copy of the full
  107,113-cell object (7.75 GB each, plus a 1.7 GB `scale.data` layer hdWGCNA never reads) to
  build a network on 7.7 % of its cells. Use `seurat_GSE243639_SNc_DA.rds` as the pipeline
  input instead — peak usage drops to ~0.66 GB across both workers. Not a hardware limit.

### Notes
- Gene count was *not* a contributor: `SetupForWGCNA` forwards `...` to `SelectNetworkGenes`,
  whose `gene_select = "variable"` default selects the 2,000 HVGs, so the TOM is ~32 MB.

---

## [0.2.0] - 2026-07-13

### Added
- `r-service/seuratdisk_compat.R` — input loading for the bridge:
  - `load_seurat()` accepts **`.rds`** (Seurat v5 native) or `.h5Seurat` (legacy), dispatching
    on file extension, and validates the path before opening it
  - `patch_seuratdisk()` — runtime shims for three SeuratDisk breakages (see *Fixed*)
- `DATA_UNZIPPED/build_seurat_GSE243639.R` — builds the GSE243639 SNc object (29 donors,
  15 PD / 14 control) from raw 10x counts: load → QC → harmony → cluster → marker-based
  cell-type labelling → `.rds`. Writes `cluster_marker_scores.csv` as an audit trail.
- `DATA_UNZIPPED/relabel_cluster7.R` — one-off correction of a mixed neuronal cluster
  mislabelled as dopaminergic.
- `.vscode/settings.json` — points the VSCode R extension at the conda R
  (`/opt/miniconda3/envs/hdWGCNA/bin/R`). **Machine-specific**; adjust when cloning.

### Fixed
- **`LoadH5Seurat` was broken for every non-spatial object.** SeuratDisk (unmaintained since
  2022) fails three ways against SeuratObject 5.x / R 4.5, all shimmed in
  `seuratdisk_compat.R`:
  1. `GetAssayData(slot=)` — made *defunct* in SeuratObject 5.0.0 (renamed `layer=`); a hard
     error, not suppressible via `lifecycle_verbosity`.
  2. `SetAssayData(slot=)` — same.
  3. `GetImages()` calls `unlist(x = assays.images, index$global$images)`; the second argument
     binds to `unlist`'s `recursive` parameter instead of being concatenated, and R ≥ 4.5
     rejects a zero-length `recursive`. Upstream paren bug, latent for years.
- **Unreadable HDF5 errors on bad input.** Passing a non-HDF5 file (an `.rds`, an `.R` script)
  dropped into the HDF5 C library and produced ~20 lines of `H5F.c` stack ending in
  `unable to read superblock`. `load_seurat()` now checks existence and extension first and
  reports one line.
- `.gitignore` did not exclude `*.rds`, `*.tar`, or `DATA_UNZIPPED/`, so a `git add -A` would
  have attempted to push ~11 GB. Build scripts under `DATA_UNZIPPED/` remain tracked.

### Changed
- `frontend/src/components/ParamForm.jsx` — file field relabelled from "h5Seurat path" to
  "Seurat path (.rds or .h5Seurat)"; placeholder now suggests `.rds`. The request parameter
  is still `h5seurat_path`, so the backend/r-service contract is unchanged.

### Notes
- **`.h5Seurat` cannot represent a Seurat v5 `Assay5`.** `SaveH5Seurat()` errors with
  "assays must have either a 'counts' or 'data' slot"; the format was abandoned before v5
  existed and no maintained writer supports it. Downgrading via
  `obj[["RNA"]] <- as(obj[["RNA"]], "Assay")` produces a loadable file, but **`.rds` is the
  supported path** — as `install.R` already noted.

---

## [0.1.1] - 2026-06-29

### Added
- Full React UI in `frontend/src/`:
  - `PipelineStepper` — horizontal 6-step indicator with break-point marker on step 4 (Soft Power)
  - `ParamForm` — parameter input form with defaults pre-filled
  - `JobStatus` — Galaxy-style status chip (running/done/failed) with 5 s polling
  - `ResultsPanel` — displays soft power PNG and kME PNG as `<img>` tags after jobs complete
  - `ModuleTable` — sortable table with colored module chips (CSS named colors from GetModules)
  - `api.js` — fetch wrappers for all 6 backend endpoints
  - `styles.css` — sidebar + results panel layout, stepper, form, table, status chip styles
- `vite.config.js` dev proxy: API paths forwarded to `localhost:8200` without CORS issues
- Break-point UX: "Run Full Pipeline" button locked until soft_power value is entered after inspecting the diagnostic plot

---

## [0.1.0] - 2026-06-29

### Added
- `r-service/plumber.R` — plumber REST bridge over hdWGCNA with endpoints:
  - `POST /test-soft-powers` — runs setup through TestSoftPowers, returns job_id
  - `POST /analyze` — full pipeline through ModuleConnectivity, requires soft_power
  - `GET /status/{job_id}` — polls in-memory job map
  - `GET /results/{job_id}/modules` — module table as JSON
  - `GET /results/{job_id}/plot` — network PNG
  - `GET /results/{job_id}/soft-power-plot` — TestSoftPowers diagnostic PNG
- `r-service/entrypoint.R` — plumber startup with `future` multisession plan
- `r-service/Dockerfile` — R 4.4 image with hdWGCNA, Seurat, WGCNA, plumber stack
- `backend/main.py` — FastAPI proxy app mirroring all r-service routes
- `backend/client.py` — thin `httpx` client for r-service communication
- `backend/Dockerfile` — Python 3.12 slim image
- `frontend/` — Vite + React 18 scaffold (pipeline stage skeleton, no component logic yet)
- `docker-compose.yml` — r-service + backend + react-ui (profile-gated) with shared volume
- `Makefile` — local dev shortcuts: `run-r`, `run-backend`, `dev`, `docker-up`, `docker-down`
- `.env.example` — environment variable template for local development
- `docs/SKILLS.md` — roadmap of planned future hdWGCNA endpoints
- `docs/UI_DESIGN.md` — wireframe and UX notes for frontend implementation
- `CLAUDE.md` — AI pair-programming context and guardrails
- `PROMPT_TEMPLATE.md` — reusable template for future work requests
