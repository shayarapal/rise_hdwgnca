# Changelog

All notable changes to this project will be documented in this file.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).  
Versioning follows [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Added
- **Soft power is now a computed number, not a guess.** `TestSoftPowers()` always computed the
  numeric scale-free-fit table, but nothing ever called `GetPowerTable()` — the only output was a
  PNG, so a power could only be chosen by eye. The pipeline now writes `soft_power_table.csv` and
  `soft_power_recommendation.json`, and `GET /results/{job_id}/soft-powers` returns the table plus
  two recommendations:
  - `recommended_power` — smallest `Power` with `SFT.R.sq >= 0.8` **and** `Power > 3`. This
    reproduces `ConstructNetwork(soft_power = NULL)`'s own rule exactly, so the recommendation and
    the network it would build can never disagree.
  - `smallest_power` — smallest `Power` with `SFT.R.sq >= 0.8`, no floor. hdWGCNA's `min_power = 3`
    means it can never auto-select powers 1–3, so this can legitimately be lower; `differ` flags it.
  On the GSE243639 DA object (`fraction = 0.05`) both resolve to **4** (R² = 0.837).
  The UI prefills the soft-power input with `recommended_power` and highlights the qualifying rows.
- `POST /gene-selection` + `GET /results/{job_id}/gene-selection[-plot]` — sweeps candidate
  `fraction` values and reports how many genes each keeps. hdWGCNA ships no tuner for `fraction`
  (there is no `TestSoftPowers` analogue), but gene selection runs on the raw object *before*
  metacells, so the whole sweep costs one object load. Each candidate calls the real
  `SetupForWGCNA` + `GetWGCNAGenes` — nothing is reimplemented — and each is isolated in a
  `tryCatch`, since `SelectNetworkGenes` errors at 0 genes and warns at ≤ 100.
- `gene_select` and `fraction` are now request parameters (were hardcoded to `"fraction"` / `0.05`).
  Both default to `None` in the backend so the r-service still owns the defaults. `"custom"` is
  rejected: it needs a `gene_list`, and this service takes no such parameter.
- `enableWGCNAThreads()` via the `WGCNA_THREADS` env var (default 4). WGCNA correlation work was
  running single-threaded. This is a **multiplier**: `entrypoint.R` runs 2 future workers, so 4
  threads is up to 8 extra R processes — set `WGCNA_THREADS=1` to disable. Falls back to
  `allowWGCNAThreads()` (thread count without a nested cluster) if the cluster cannot be created.
- `pipeline_libs()` now attaches `WGCNA`, `tidyverse`, `cowplot`, `igraph`, `enrichR` to match the
  reference analysis preamble. `enrichR` is attached only; no enrichment is run.
- The soft-power plot is also written as `soft_power_plot.pdf` (`ggsave`, 10×5), reachable through
  the file browser. The plot now marks the recommended power via `PlotSoftPowers(selected_power=)`.
- **`DATA_GSE233866/build_vta_intact_network.R` + `run_dme_vta.R`** — the VTA counterpart to the
  existing SNc `intact`/`lesioned`/combined-DME trio, which previously had no equivalent (VTA only
  had a standalone `lesioned` network, built solely as input to the SNc-vs-VTA DME test). Finds
  CACNA1D's module responds to 6-OHDA lesioning in VTA in the **same direction** as SNc (module
  eigengene down, kME up: 0.230→0.283) but at far lower magnitude (p.adj=2.2×10⁻⁵ vs. SNc's
  6.1×10⁻²⁶⁹ — roughly 264 orders of magnitude less extreme). This means "only SNc responds to
  lesioning" is **not** supported by this project's data and has been removed from every doc that
  implied it; the region-selectivity claim now rests on the direct SNc-vs-VTA-during-lesioning
  comparison (`snc_vs_vta_combined_dme_lesioned/`, unaffected, p.adj=4.5×10⁻¹²³), not on VTA
  showing no response at all. Full writeup: `analysis_results/mouse_GSE233866/README.md` §4.
- `archive/human_GSE243639/` — the built human Seurat objects and analysis results were archived
  (not deleted) to allow a clean rebuild from raw data; `DATA_UNZIPPED/` still has the raw counts
  and every build script needed to reproduce them.
- Per-folder `STATS.md` + `INTERPRETATION.tex` next to every `analysis_results/mouse_GSE233866/`
  result folder (10 folders total), sourced directly from each folder's own CSV/JSON rather than
  copied from a central summary — see the "Removed" note below for why that matters here.
- Published an interactive "hdWGCNA Field Guide" reference explaining every statistic type
  (soft power, kME, module eigengenes, preservation Zsummary, DME log2FC/p.adj, Enrichr output)
  produced by this pipeline, using the mouse SNc network as a worked example.

### Fixed
- **`ConstructNetwork(soft_power = NULL)` silently yields `Inf`** when no power clears the
  scale-free threshold: its rule is `subset(...) %>% .$Power %>% min`, and `min(numeric(0))` is
  `Inf`, which it passes into WGCNA without complaint. `recommend_power()` returns `null` plus an
  explicit warning instead, and `/analyze` + `/module-preservation` now reject a non-finite
  `soft_power` at the boundary.
- `/results/{job_id}/soft-powers` uses `@serializer json list(na = "null")` rather than
  `unboxedJSON`, whose `toJSON()` defaults encode `NA_integer_` as the **string `"NA"`** and `NULL`
  as `{}` — either of which would reach the client as a truthy value in exactly the case
  (no qualifying power) the caller most needs to detect.
- **`r-service/Dockerfile` never copied `seuratdisk_compat.R`**, which `plumber.R` sources at
  startup — the container could not boot.
- **`frontend/vite.config.js` did not proxy `/module-preservation`**, so comparison mode was broken
  under `make run-ui`. `frontend/nginx.conf` was missing `module-preservation` and `files` too.
- `ConstructNetwork` now passes `overwrite_tom = TRUE`. It defaulted to `FALSE`, reusing any TOM
  already in `out_dir`; now that `fraction` can change the gene set, a stale TOM is not merely old
  but built on **different genes**, which would silently produce modules for genes nobody asked for.
- `SetupForWGCNA` errors if gene selection yields 0 genes, naming the parameter responsible.
  `gene_select`/`fraction` are not `SetupForWGCNA` formals — they ride through `...` to
  `SelectNetworkGenes` — so this also trips if a future hdWGCNA stops forwarding them.
- **All five mouse GSE233866 network build scripts queried the wrong species' KEGG library** —
  `dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")`, copy-pasted unchanged from the
  human GSE243639 scripts into every mouse one (`build_snc_network.R`, `build_vta_network.R`,
  `build_snc_intact_network.R`, `build_snc_lesioned_network.R`, `build_vta_lesioned_network.R`).
  Switched to `KEGG_2019_Mouse` and reran all five networks against freshly re-downloaded raw
  counts; every non-enrichment output (`modules.csv`, soft power) was verified byte-identical to
  the pre-fix run, confirming only the KEGG enrichment results were affected. The previously-cited
  SNc "Dopaminergic synapse" hit holds up under the correct library at essentially the same p.adj
  (2.8×10⁻⁴ → 2.77×10⁻⁴).
- **CACNA1D's healthy-baseline kME was misread from the wrong `kME_*` column** in both
  `healthy_baseline/{snc_network,vta_network}/modules.csv` — reported as 0.048 (SNc, "very weak")
  and 0.094 (VTA, "weak, stronger than SNc"); actually 0.442 (SNc) and 0.296 (VTA), both moderate,
  with SNc — not VTA — the stronger connection, the reverse of what was originally stated. This
  had already propagated into the mouse README, `SYNTHESIS.md`'s cross-study table, and the
  published Field Guide before being caught; all were corrected with an inline note showing the
  old value, the new value, and why, rather than a silent rewrite. The human GSE243639 kME figures
  were independently re-checked against their own source files and are correct — this misread was
  isolated to the two mouse healthy-baseline files.

### Removed
- `analysis_results/mouse_GSE233866/healthy_baseline/STATISTICS_SUMMARY.md` — a hand-maintained
  combined-numbers table that had already drifted out of sync with the per-folder `STATS.md` files
  and `README.md` (missing study #8 entirely) days after being written. Removed rather than kept
  in sync going forward; the per-folder `STATS.md` files sit next to their source data and
  `README.md` is the single narrative summary.

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
  the data, not a bug. **`/module-preservation` now writes `donor_counts.csv` on every run so
  this n sits next to the result rather than being buried.**

---

## [0.5.0] - 2026-07-13

> **The bridge can now answer the question the project was built for.** `/analyze` only ever
> produced a single pooled network; comparing PD against control required
> `ModulePreservation`, which was unwrapped. It now is.

### Added
- **`POST /module-preservation`** — reference/query module preservation, wrapping hdWGCNA's
  `ProjectModules` + `ModulePreservation`. Splits the input object on `condition_col`, builds the
  full network in `ref_group`, projects those modules into `query_group`, and scores whether each
  one survives. Same pattern as the existing endpoints: file-path in, `job_id` out, poll, fetch.
  - Required: the six base fields plus `soft_power`, `condition_col`, `ref_group`, `query_group`
  - Optional: `n_permutations` (default **250**, the value used in the hdWGCNA tutorial — the
    `ModulePreservation()` default is 500), `preservation_name` (default `"<ref>-vs-<query>"`)
  - Call order follows the [tutorial](https://smorabit.github.io/hdWGCNA/articles/module_preservation.html)
    exactly: `ProjectModules` runs **before** the query's metacells, and `SetDatExpr` runs on the
    reference **before** the query. `ProjectModules` invokes `SetupForWGCNA` and `ModuleEigengenes`
    on the query internally, so neither is called on it directly.
- `GET /results/{job_id}/preservation` — per-module preservation statistics as JSON, sorted
  **least-preserved first**. Zsummary convention (Langfelder et al. 2011): `> 10` strongly
  preserved, `2–10` weak, `< 2` not preserved.
- `GET /results/{job_id}/preservation-plot` — `PlotModulePreservation(statistics = "summary")` PNG.
- `GET /results/{job_id}/donor-counts` — cells per donor per condition, and which clear
  `MetacellsByGroups(min_cells = 100)`. Written by `/module-preservation`; see *Known issues*.
- `condition_col` + `ref_group` are now **optional** on `POST /test-soft-powers`. When given, it
  subsets to the reference condition first, so the soft power is chosen from the curve of the
  network `/module-preservation` will actually build — not from the pooled object.
- `frontend/` — "Compare conditions" toggle in `ParamForm`, `PreservationTable` (Zsummary colour-
  banded by the convention above), preservation plot + table in `ResultsPanel`. In comparison mode
  the second-phase button calls `/module-preservation` instead of `/analyze`.

### Changed
- **De-duplicated the pipeline block in `plumber.R`.** `run_setup_through_soft_powers()` and
  `run_full_pipeline()` held the same ~40 lines twice, which is why both the `gene_select`/`fraction`
  fix and the `set.seed(42)` fix in 0.4.0 had to be applied in two places. Now
  `setup_and_metacells()` / `test_soft_powers()` / `build_network()` / `write_module_outputs()`,
  shared by all three job entrypoints. The result endpoints likewise share `serve_csv()` /
  `serve_png()`. No hdWGCNA call, argument, or default changed.
- `/module-preservation` reuses the existing result filenames (`modules.csv`, `network_plot.png`,
  `soft_power_plot.png`) for the **reference** network, so `/results/{job_id}/modules`, `/plot` and
  `/soft-power-plot` work against a preservation job with no new code.

### Fixed
- **`GET /results/{job_id}/modules` was double-encoding its response.** The handler called
  `toJSON(df, ...)` *and* declared `@serializer unboxedJSON`, so plumber JSON-encoded the
  already-JSON string: clients received a quoted string (`"[{\"gene\":...}]"`) instead of an array
  of row objects. `frontend/src/App.jsx` guards with `Array.isArray(rows) ? rows : []`, so the
  module table silently rendered **empty** rather than erroring. The handler now returns the
  data.frame and lets the serializer encode it once. Caught by calling the live endpoint; the
  phase-2 tests mock the R service and so could not see it.
- README endpoint table said the backend listens on **8080**; `docker-compose.yml`, the `Makefile`
  and the Vite proxy all say **8200**.
- A missing metadata column surfaced as `all arguments must have the same length` from deep inside
  `as.data.frame(table(...))`. `run_module_preservation()` now validates `condition_col`,
  `cell_type_col` and `group_by` up front and names the offending column plus the ones available;
  a bad `ref_group`/`query_group` likewise now lists the values actually present.

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
