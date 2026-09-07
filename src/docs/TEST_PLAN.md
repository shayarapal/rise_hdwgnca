# Test Plan — hdWGCNA Bridge Service

Three phases, each a superset of the previous in terms of infrastructure required.
Run them in order; a failing Phase 1 is a blocker for Phase 2.

---

## Phase 1 — Static Analysis (no running services)

**Goal:** Prove the repo is internally consistent without building or running anything.

Tests (`tests/phase1/test_static.py`):

| Test | What it checks |
|---|---|
| `test_critical_files_exist` | All required source files are present on disk |
| `test_dockerfile_port_matches_entrypoint` | `EXPOSE` in `r-service/Dockerfile` matches the default PORT in `entrypoint.R` |
| `test_plumber_endpoints_defined` | `plumber.R` contains all six expected `@post`/`@get` decorators |
| `test_docker_compose_topology` | Services named `r-service`, `backend`, `react-ui`; correct host ports; shared volume wired to both |
| `test_r_service_url_env_var` | `backend/client.py` reads `R_SERVICE_URL`; `docker-compose.yml` sets it to `http://r-service:8100` |
| `test_python_syntax` | `py_compile` passes on every `.py` file in `backend/` |
| `test_pipeline_request_required_fields` | Pydantic raises `ValidationError` when any required field is omitted |
| `test_analyze_soft_power_guard` | `soft_power=None` triggers a 422 in the FastAPI `/analyze` route (unit-level, no HTTP server) |

**Run:** `make test-phase1`
**Infrastructure:** none — pure Python + file I/O.

---

## Phase 2 — Backend HTTP Smoke Tests (FastAPI + mocked R service)

**Goal:** Verify the FastAPI layer enforces its contract and correctly proxies to the R service,
without needing Docker or a real R process. Uses `httpx.MockTransport` to fake the R service.

Tests (`tests/phase2/test_backend.py`):

| Test | What it checks |
|---|---|
| `test_post_test_soft_powers_missing_fields` | Backend returns 400 when required fields are absent |
| `test_post_analyze_missing_soft_power` | Backend returns 422 when `soft_power` is omitted |
| `test_post_analyze_proxies_body` | Body forwarded to R service verbatim (including `soft_power`) |
| `test_get_status_proxies` | `/status/{job_id}` passes job_id through and returns R service response |
| `test_get_modules_proxies` | `/results/{job_id}/modules` returns JSON from R service |
| `test_get_plot_content_type` | `/results/{job_id}/plot` returns `image/png` |
| `test_get_soft_power_plot_content_type` | `/results/{job_id}/soft-power-plot` returns `image/png` |
| `test_r_service_error_propagates` | Non-2xx from R service bubbles up as HTTPException with matching status code |

**Run:** `make test-phase2`
**Infrastructure:** Python only — `pytest` + `httpx` (already in `requirements.txt`) + `fastapi[testclient]`.

---

## Phase 3 — Docker Compose Integration Tests (full stack, no real data)

**Goal:** Prove the two containers start, talk to each other, and enforce their API contracts
end-to-end. Does **not** require a real `.h5Seurat` file — all pipeline calls are
expected to fail with a file-not-found error from inside R, which is the correct
observable behaviour given a fake path.

Tests (`tests/phase3/test_integration.py`):

| Test | What it checks |
|---|---|
| `test_r_service_reachable` | `GET http://localhost:8100/` returns *something* (even a 404 is fine) |
| `test_backend_reachable` | `GET http://localhost:8200/` returns *something* |
| `test_backend_missing_fields_400` | `POST /test-soft-powers` with empty body → 422 from backend |
| `test_backend_analyze_missing_soft_power_422` | `POST /analyze` without `soft_power` → 422 |
| `test_r_service_missing_fields_400` | Direct `POST http://localhost:8100/test-soft-powers` with empty body → 400 |
| `test_unknown_job_status_404` | `GET /status/does-not-exist` → 404 |
| `test_fake_path_job_lifecycle` | Submit job with fake h5seurat path → get `job_id` → poll until `failed` → error message non-empty |

**Run:** `make test-phase3`  
**Infrastructure:** `docker compose up --build` must be running before executing this phase.
The Makefile target does **not** start the stack — run `make docker-up` first, wait for
both services to be healthy, then run `make test-phase3`.

---

## Open Questions

Before Phase 3 can produce a full-pipeline green test, you need to supply:

- A real `.h5Seurat` file accessible under `/shared/data/` in the container
- The exact metadata column name for cell type (`cell_type_col`)
- The `group_by` list and `group_name` value that match your data
- A confirmed `soft_power` integer from a previous `TestSoftPowers` run

These are inputs, not assumptions — no placeholder values will be added to the test suite.
