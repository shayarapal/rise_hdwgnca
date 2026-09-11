# hdWGCNA Bridge

A minimal `plumber`-based R microservice that exposes the [hdWGCNA](https://github.com/smorabit/hdWGCNA)
single-cell co-expression pipeline over HTTP, so it can be called from a Python
backend exactly the way you'd call a native Python library.

## Why this exists

[PyWGCNA](https://github.com/mortazavilab/PyWGCNA) gives bulk RNA-seq WGCNA users a
native Python interface. **hdWGCNA has no Python equivalent** — it is R-only and
requires Seurat objects. This bridge does not reimplement hdWGCNA; it wraps the
existing R functions in HTTP endpoints so a Python application can drive the same
pipeline without running R itself.

```
Python backend  ──HTTP──▶  R plumber bridge  ──calls──▶  hdWGCNA / Seurat (R)
```

## Scope

This service wraps the following hdWGCNA pipeline stages, in the order used in the
[official tutorial](https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html):

| Stage | hdWGCNA function |
|---|---|
| Setup | `SetupForWGCNA` |
| Metacells | `MetacellsByGroups`, `NormalizeMetacells` |
| Expression matrix | `SetDatExpr` |
| Soft power selection | `TestSoftPowers` |
| Network construction | `ConstructNetwork` |
| Module detection | `ModuleEigengenes`, `ModuleConnectivity` |
| Condition comparison | `ProjectModules`, `ModulePreservation` |

Nothing beyond this list is in scope unless explicitly added.

## Endpoints

### r-service (port 8100)

| Method | Path | Description |
|---|---|---|
| `POST` | `/test-soft-powers` | Runs setup → metacells → SetDatExpr → TestSoftPowers; returns job_id |
| `POST` | `/analyze` | Runs full pipeline; requires `soft_power` in request body |
| `POST` | `/module-preservation` | Builds the network in `ref_group`, projects it into `query_group`, scores preservation |
| `POST` | `/gene-selection` | Sweeps candidate `fraction` values, reporting genes kept by each; returns job_id |
| `GET` | `/status/{job_id}` | Returns job status (`running` / `done` / `failed`) |
| `GET` | `/results/{job_id}/modules` | Returns the module table as JSON |
| `GET` | `/results/{job_id}/soft-powers` | Numeric scale-free-fit table + the recommended soft power |
| `GET` | `/results/{job_id}/gene-selection` | Genes kept per candidate fraction, as JSON |
| `GET` | `/results/{job_id}/gene-selection-plot` | Genes-vs-fraction curve (PNG) |
| `GET` | `/results/{job_id}/plot` | Returns the module network plot (PNG) |
| `GET` | `/results/{job_id}/soft-power-plot` | Returns the TestSoftPowers diagnostic plot (PNG) |
| `GET` | `/results/{job_id}/preservation` | Per-module preservation statistics as JSON, least-preserved first |
| `GET` | `/results/{job_id}/preservation-plot` | Returns the PlotModulePreservation summary plot (PNG) |
| `GET` | `/results/{job_id}/donor-counts` | Cells per donor per condition, and which clear `min_cells = 100` |

For a `/module-preservation` job, `/modules`, `/plot` and `/soft-power-plot` describe the
**reference** network.

### backend (port 8200)

Mirrors all r-service endpoints exactly. Acts as a stable API surface for the
frontend and any other consumers.

## Typical workflow

```
0. POST /gene-selection    →  { job_id }   (optional: tune `fraction` before anything else)
   GET  /results/{job_id}/gene-selection   →  genes kept per candidate fraction

1. POST /test-soft-powers  →  { job_id }
2. GET  /status/{job_id}   →  poll until "done"
3. GET  /results/{job_id}/soft-powers      →  numeric table + recommended_power
4. POST /analyze  { ..., soft_power: 4 }   →  { job_id }
5. GET  /status/{job_id}   →  poll until "done"
6. GET  /results/{job_id}/modules          →  JSON module table
7. GET  /results/{job_id}/plot             →  network PNG
```

## Choosing the soft power

`GET /results/{job_id}/soft-powers` returns the numeric scale-free-fit table (the same one
`PlotSoftPowers` draws) plus two recommendations, so the power can be chosen from data rather
than by eye:

```json
{
  "table": [{ "Power": 4, "SFT.R.sq": 0.837, "slope": -17.99, "mean.k.": 495.2 }],
  "recommended_power": 4,
  "smallest_power": 4,
  "differ": false,
  "sft_threshold": 0.8,
  "min_power": 3,
  "max_sft_r_sq": 0.9895,
  "warning": null
}
```

- **`recommended_power`** — the smallest power with `SFT.R.sq >= 0.8` **and** `Power > 3`. This
  reproduces `ConstructNetwork(soft_power = NULL)`'s own rule exactly, so the number you are shown
  is the network hdWGCNA would build. This is what the UI prefills.
- **`smallest_power`** — the smallest power with `SFT.R.sq >= 0.8`, with no floor. hdWGCNA's
  `min_power = 3` means it can never auto-select powers 1–3, so this is occasionally lower;
  `differ` is `true` when the two disagree.
- Both are **`null`** — never `Inf` — when no power reaches the threshold, and `warning` says so.
  (Left to itself, `ConstructNetwork(soft_power = NULL)` computes `min(numeric(0))` = `Inf` here
  and passes it into WGCNA silently.)

## Choosing the gene fraction

`gene_select` / `fraction` control which genes enter the network (default `"fraction"` / `0.05`,
the tutorial values — a gene is kept if it is detected in at least 5% of cells; this is *not* the
top 5% of genes). hdWGCNA provides no tuner for this, so `POST /gene-selection` sweeps candidate
fractions and reports how many genes each keeps, letting you pick a gene-set size deliberately:

| fraction | genes kept (GSE243639 DA neurons, 28,437 total) |
|---|---|
| 0.01 | 11,958 |
| 0.05 | 7,465 ← default |
| 0.10 | 5,011 |
| 0.20 | 2,518 |

Too high and real co-expression partners are dropped; too low and the network fills with noise.
hdWGCNA warns below 100 genes and errors at 0.

## Request body

All POST endpoints accept the same JSON shape (`soft_power` is only required for `/analyze`
and `/module-preservation`):

```json
{
  "h5seurat_path": "/shared/data/brain.h5Seurat",
  "out_dir":       "/shared/results/run1",
  "cell_type_col": "cell_type",
  "group_by":      ["cell_type", "Sample"],
  "group_name":    "Excitatory Neurons",
  "wgcna_name":    "tutorial",
  "k":             25,
  "max_shared":    15,
  "network_type":  "signed",
  "soft_power":    9
}
```

`group_by` must be exactly two columns, and the **second** is the batch variable — it is passed
to `ModuleEigengenes(group.by.vars = ...)` and used as the donor column in `donor_counts.csv`.

## Comparing two conditions

`/module-preservation` answers "does this module still hold together in the other condition?"
It takes the fields above plus:

```json
{
  "condition_col": "condition",
  "ref_group":     "control",
  "query_group":   "PD",
  "n_permutations": 250
}
```

The network is built in `ref_group` and tested in `query_group`. `n_permutations` defaults to
250 (the value used in the hdWGCNA tutorial; the `ModulePreservation()` default is 500) and is
the dominant cost — this endpoint is far slower than `/analyze`.

Pass `condition_col` + `ref_group` to `/test-soft-powers` as well, so the soft power you pick
comes from the reference condition's curve rather than the pooled object's.

Reading the result: find your gene of interest in `/results/{job_id}/modules` to get its module,
then look that module up in `/results/{job_id}/preservation`. By the usual convention
(Langfelder et al. 2011), `Zsummary > 10` means strongly preserved, `2–10` weak, and `< 2` not
preserved — a **low** Zsummary means the module's co-expression structure breaks down in the
query condition. Check `/results/{job_id}/donor-counts` for the n behind that number.

## Requirements

- R 4.4+
- `plumber`, `future`, `promises`, `uuid`, `jsonlite`
- `Seurat`, `SeuratDisk`, `WGCNA`, `igraph`, `tidyverse`, `ggraph`, `patchwork`
- `hdWGCNA` (installed from `smorabit/hdWGCNA`, `dev` branch)
- Bioconductor: `UCell`, `GenomicRanges`, `GeneOverlap`

See [r-service/Dockerfile](r-service/Dockerfile) for the exact install sequence.

> **Note on SeuratDisk**: `SeuratDisk` (used for `LoadH5Seurat`) has not been
> actively maintained since 2022. If your `.h5Seurat` files were written by a
> recent Seurat v5 build, you may need to switch to `readRDS` + `.rds` files.
> See the comment in `src/r-service/Dockerfile`.

## Running

### Docker (recommended)

```bash
docker compose up --build r-service backend
```

Services listen on `:8100` (r-service) and `:8200` (backend).

### Locally (without Docker)

```bash
cp .env.example .env          # edit SHARED_DIR to point at your data (outside this repo)
make setup-local              # creates $SHARED_DIR/{DATA_GSE233866,DATA_UNZIPPED,results}
make install                  # installs R packages, pip deps, node_modules (once, takes ~1-3 h for R)

make run-r          # starts plumber on :8100  (terminal 1)
make run-backend    # starts FastAPI on :8200  (terminal 2)
make run-ui         # starts Vite dev on :4000 (terminal 3)
# or
make dev            # starts r-service + backend in parallel; open a third terminal for run-ui
```

### Quick smoke test

```bash
curl -X POST http://localhost:8100/test-soft-powers \
  -H "Content-Type: application/json" \
  -d '{
    "h5seurat_path": "/shared/data/brain.h5Seurat",
    "out_dir":       "/shared/results/run1",
    "cell_type_col": "cell_type",
    "group_by":      ["cell_type", "Sample"],
    "group_name":    "Excitatory Neurons",
    "wgcna_name":    "tutorial",
    "k":             25,
    "max_shared":    15,
    "network_type":  "signed"
  }'
# → {"job_id":"<uuid>"}

curl http://localhost:8100/status/<uuid>
# → {"status":"running"} ... {"status":"done"}

curl http://localhost:8100/results/<uuid>/soft-power-plot --output sp.png
```

## What this is not

- Not a Python port of hdWGCNA. There is no Python implementation of the algorithm
  anywhere in this repo — only an HTTP wrapper around the R package.
- Not a general-purpose R API. Endpoints are added only for pipeline stages
  actually needed by the consuming application.
- Not a replacement for reading the
  [hdWGCNA documentation](https://smorabit.github.io/hdWGCNA/) — parameter meaning,
  defaults, and biological interpretation come from there, not from this service.

## Use of AI assistance

The **software** in this repository was written with AI assistance. Claude Code was used
as a pair-programming tool to build the plumber bridge, the FastAPI proxy, the Docker
setup, and the test suite, working from the official hdWGCNA tutorials. `CLAUDE.md`
records the guardrails it was held to — chiefly that hdWGCNA functions are called exactly
as documented upstream, and that no part of the algorithm is reimplemented,
approximated, or ported to Python.

The **scientific results** were not generated by AI. Every number in the analysis outputs
— module assignments, kME values, soft-power selections, differential module eigengene
statistics, preservation Z-scores, and enrichment tables — is computed by hdWGCNA, WGCNA,
Seurat and Enrichr running on public GEO data (GSE233866). A language model does not
produce these values and cannot: they come out of the R packages themselves, and rerunning
the pipeline from the raw counts reproduces them. The AI's role stopped at writing the
code that calls those packages and at helping check the resulting numbers back against the
output files.

Put simply: AI helped build the instrument; the measurements are the instrument's.

## Repo layout

Everything except the analysis outputs lives under `src/`. The analysis outputs themselves
live outside the repository, alongside the datasets — see "Where the data lives" below.

```
src/
  r-service/          # plumber bridge
    plumber.R         # endpoint definitions
    entrypoint.R      # startup: loads plumber, sets future plan
    Dockerfile
  backend/            # Python FastAPI proxy calling the bridge
    main.py
    client.py
    Dockerfile
    requirements.txt
  frontend/           # React UI scaffold (component code not yet implemented)
    src/App.jsx       # pipeline stage skeleton
    package.json
  tests/              # phase1 (static), phase2 (mocked), phase3 (integration)
  docs/
    SKILLS.md         # planned future endpoint capabilities
    UI_DESIGN.md      # wireframe and UX notes for frontend implementation
    CHANGELOG.md
    AUDIT_FINDINGS.md
    PROMPT_TEMPLATE.md            # reusable template for future work requests
    HDWGCNA_SYSTEM_CHECK_PROMPT.md

reminders/            # working notes

docker-compose.yml    # service topology (build contexts point into src/)
Makefile              # local dev shortcuts
.env.example
CLAUDE.md             # AI pair-programming context and guardrails
```

### Where the data and the analysis scripts live

**This repository is the hdWGCNA HTTP bridge only.** The datasets *and* the per-study R
analysis scripts both live **outside** it, in the directory named by `SHARED_DIR` (see
`.env.example`), which is bind-mounted at `/shared` in both containers:

```
$SHARED_DIR/                     # e.g. C:/projects/rise_hdwgnca_data — NOT in the repo
  DATA_GSE233866/                # mouse: counts, built Seurat objects, run outputs, *.R
  DATA_UNZIPPED/                 # human: GSE243639_RAW.tar, counts, SEURAT_INPUT/, *.R
  analysis_results/              # pipeline outputs: modules, kME, DME, enrichment, figures
  results/                       # earlier exploratory runs
```

`analysis_results/` holds every figure and table behind the write-up — per-network `STATS.md`,
`SYNTHESIS.md`, module CSVs and Enrichr exports. It was moved out of the repository on
2026-09-06 and is no longer tracked; earlier revisions remain in git history
(`git checkout 446ce5f -- analysis_results/`). No code reads it — it is output only.

The 24 R scripts that produced the published analysis (Seurat build, network construction,
DME, module preservation, diagnostics) sit next to the data they read, in those same two
subdirectories. They were untracked here on 2026-09-07; earlier revisions remain in git
history: `git checkout 4e8e769 -- DATA_GSE233866/ DATA_UNZIPPED/`.

Those two subdirectory names are **load-bearing**: the scripts address data by absolute
container path, e.g. `/shared/DATA_GSE233866/seurat_GSE233866_SNc_VTA.rds`. Renaming them
breaks every script.

To set this up from a fresh clone:

```bash
cp .env.example .env          # then edit SHARED_DIR to a path outside the repo
make setup-local              # creates the subdirectory skeleton
# download raw counts from GEO: GSE233866 (mouse), GSE243639 (human)
# place the analysis .R scripts in the matching subdirectories
```

## Code and data availability

### Data

This project generated no new data. All analyses use publicly available datasets
downloaded from the NCBI Gene Expression Omnibus:

| Accession | Organism | Content |
|---|---|---|
| `GEO: GSE233866` | Mouse | SNc/VTA dopaminergic atlas; untreated baseline and 6-OHDA lesion/intact counts |
| `GEO: GSE243639` | Human | SNc dopaminergic neurons, PD and control donors |

Neither dataset was produced by this project, and neither is redistributed here — see
"Where the data and the analysis scripts live" above for the expected local layout.

### Code

The bridge service in this repository is released under GPL-3.0-or-later (see `LICENSE`),
matching the license of hdWGCNA, which it loads in-process.

The per-study R analysis scripts that produced the published results are **not** tracked in
this repository. They live alongside the data in `$SHARED_DIR`, and earlier revisions remain
in git history at `git checkout 4e8e769 -- DATA_GSE233866/ DATA_UNZIPPED/`.

### Environment used for the published analysis

Results were produced with the `r-service` image built 2026-07-29. The GitHub dependencies are
pinned to the exact commits in that image, so rebuilding reproduces the same stack:

| Component | Version | Pin |
|---|---|---|
| R | 4.4.3 | `rocker/r-ver:4.4.3` |
| hdWGCNA | 0.4.12 | `e3344d1f7bbac4264adf94f7aa31e0802fa8282d` |
| WGCNA | 1.74 | CRAN (RSPM snapshot fixed by the base image) |
| Seurat | 5.5.1 | CRAN (as above) |
| SeuratObject | 5.4.0 | CRAN (as above) |
| SeuratDisk | 0.0.0.9021 | `877d4e18ab38c686f5db54f8cd290274ccdbe295` |
| igraph | 2.3.3 | CRAN (as above) |
| enrichR | 3.4 | CRAN (as above) |
| ggforestplot | 0.1.0 | `547617e63fa481a5f28ffc56c07d46be4af688b2` |

Do not float `hdWGCNA` back to `ref='dev'`. Module assignments are version-sensitive, and the
published numbers came from the pinned commit.

Note that `src/r-service/install.R` (the non-Docker path) installs CRAN packages from
`cloud.r-project.org` and therefore resolves to current versions rather than the snapshot
above. Use the Docker path to reproduce the published environment.
