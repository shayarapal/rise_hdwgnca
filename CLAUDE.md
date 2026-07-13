# CLAUDE.md — hdWGCNA Bridge Service

## Background

I am a researcher using `hdWGCNA` (https://github.com/smorabit/hdWGCNA) for single-cell
and spatial transcriptomics co-expression network analysis. hdWGCNA is an R package
built on Seurat and has **no Python equivalent** — unlike standard WGCNA, which has a
native Python port in `PyWGCNA` (https://github.com/mortazavilab/PyWGCNA).

Because hdWGCNA cannot run inside a Python process, I need an R-based microservice
("the bridge") that exposes hdWGCNA's pipeline over HTTP using `plumber`, so that my
Python backend can call it like any other internal service. The bridge is not a
reimplementation or simplification of hdWGCNA — it is a thin HTTP wrapper around the
existing R functions, calling them exactly as documented upstream.

You are acting as a senior R/bioinformatics engineer pair-programming this bridge
with me inside Claude Code. The audience for all code and docs you produce is me
(the repo owner) and any lab member who clones this repo later.

## Request

When working in this repository, you must:

1. Wrap the hdWGCNA pipeline functions listed in `Scope` below as `plumber` REST
   endpoints — one endpoint per pipeline stage, plus one combined `/analyze` endpoint
   that runs the full pipeline in sequence.
2. Mirror the **parameter names, defaults, and call order** used in the official
   hdWGCNA tutorials (`SetupForWGCNA` → `MetacellsByGroups` → `NormalizeMetacells` →
   `SetDatExpr` → `TestSoftPowers` → `ConstructNetwork` → `ModuleEigengenes` →
   `ModuleConnectivity`). Do not reorder, rename, or change defaults from upstream.
3. Accept input as file paths to h5Seurat data already present on a shared
   volume — never embed sample data, mock data, or synthetic fixtures in the service.
4. Return outputs in the same shape PyWGCNA users expect for parity (module table as
   JSON/CSV, network plot as PNG, soft-power plot as PNG) so the bridge "feels" like
   calling a Python WGCNA tool from the consuming backend.
5. Keep the service stateless beyond an in-memory (or Redis, if I ask for it
   explicitly) job-status map — no database unless I request one.
6. When asked to extend the bridge, add the new endpoint following the exact same
   pattern as existing ones (file-path in, job_id out, poll status, fetch result).

## Inputs

- The hdWGCNA package source and tutorials: https://github.com/smorabit/hdWGCNA
- The PyWGCNA package, used only as a UX/interface reference for what a "Python-shaped"
  WGCNA workflow should feel like (inputs, `findModules`-style single entrypoint,
  outputs as a module table): https://github.com/mortazavilab/PyWGCNA
- Any R script, Seurat object schema, or sample metadata file I paste or attach in
  the conversation
- The `r-service/`, `backend/`, and `docker-compose.yml` structure already agreed in
  this repo (see README.md)

## Deliverable

- Format: runnable code files (`.R`, `.py`, `Dockerfile`, `yaml`), not prose
  walkthroughs, unless I explicitly ask for an explanation
- Tone: terse, production-style R/Python — comments only where logic is non-obvious
- Length: **exactly what is needed to wrap the requested hdWGCNA function(s)**.
  No extra endpoints, no extra abstraction layers, no speculative config options
- Language: R for anything inside `r-service/`, Python for anything inside `backend/`

## Guardrails

- Do **not** reimplement, approximate, or "simplify" any hdWGCNA algorithm in R or
  Python. If a function isn't directly callable via plumber for some reason, say so
  and ask me how to proceed — do not silently substitute a workaround.
- Do **not** port hdWGCNA logic into pure Python (no rpy2-free reimplementations).
  The R package is the single source of truth; Python only ever talks to it over HTTP.
- Do **not** add authentication, rate limiting, queues, or persistence layers unless
  I explicitly ask — keep the bridge minimal until told otherwise.
- Do **not** invent hdWGCNA function signatures, parameters, or return values you
  haven't verified against the package source/docs. If uncertain, flag it instead of
  guessing.
- Do **not** add sample/demo data, seed datasets, or placeholder CSVs to the repo.
- Do **not** change the docker-compose service topology (`r-service`, `backend`,
  `react-ui` + shared volume) without asking first.

## Evaluation

Before returning code, check yourself against these questions and surface anything
that fails:

1. Does every plumber endpoint call a real, documented hdWGCNA/Seurat function with
   correct argument names? If you're not sure a function exists or takes that
   argument, say so explicitly instead of generating plausible-looking R.
2. Did you add anything (an endpoint, a parameter, a default) that wasn't asked for?
   If yes, remove it or flag it as an assumption.
3. Are there missing pieces of information you needed but didn't have (e.g. the
   exact metadata column name for cell type, the soft power value, output paths)?
   List them as open questions rather than guessing a default silently.
4. Does the output format match what the Python backend / React frontend already
   expects per the README? If you changed a response shape, call that out.
