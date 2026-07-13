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

Nothing beyond this list is in scope unless explicitly added.

## Endpoints

### r-service (port 8100)

| Method | Path | Description |
|---|---|---|
| `POST` | `/test-soft-powers` | Runs setup → metacells → SetDatExpr → TestSoftPowers; returns job_id |
| `POST` | `/analyze` | Runs full pipeline; requires `soft_power` in request body |
| `GET` | `/status/{job_id}` | Returns job status (`running` / `done` / `failed`) |
| `GET` | `/results/{job_id}/modules` | Returns the module table as JSON |
| `GET` | `/results/{job_id}/plot` | Returns the module network plot (PNG) |
| `GET` | `/results/{job_id}/soft-power-plot` | Returns the TestSoftPowers diagnostic plot (PNG) |

### backend (port 8080)

Mirrors all r-service endpoints exactly. Acts as a stable API surface for the
frontend and any other consumers.

## Typical workflow

```
1. POST /test-soft-powers  →  { job_id }
2. GET  /status/{job_id}   →  poll until "done"
3. GET  /results/{job_id}/soft-power-plot  →  inspect PNG, choose soft_power value
4. POST /analyze  { ..., soft_power: 9 }  →  { job_id }
5. GET  /status/{job_id}   →  poll until "done"
6. GET  /results/{job_id}/modules          →  JSON module table
7. GET  /results/{job_id}/plot             →  network PNG
```

## Request body

Both POST endpoints accept the same JSON shape (`soft_power` is only required for `/analyze`):

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
> See the comment in `r-service/Dockerfile`.

## Running

### Docker (recommended)

```bash
docker compose up --build r-service backend
```

Services listen on `:8100` (r-service) and `:8200` (backend).

### Locally (without Docker)

```bash
cp .env.example .env          # edit SHARED_DIR to point at your data
make setup-local              # creates local_data/data and local_data/results
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

## Repo layout

```
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
docs/
  SKILLS.md         # planned future endpoint capabilities
  UI_DESIGN.md      # wireframe and UX notes for frontend implementation
docker-compose.yml
Makefile            # local dev shortcuts
.env.example
CHANGELOG.md
CLAUDE.md           # AI pair-programming context and guardrails
PROMPT_TEMPLATE.md  # reusable template for future work requests
```
