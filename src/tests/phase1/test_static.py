"""
Phase 1 — Static Analysis
No running services required. Pure file I/O and import-level checks.
"""
import pathlib
import py_compile
import re
import sys

import pytest
import yaml

ROOT = pathlib.Path(__file__).parent.parent.parent


# ---------------------------------------------------------------------------
# 1. Required files exist
# ---------------------------------------------------------------------------

REQUIRED_FILES = [
    "r-service/plumber.R",
    "r-service/entrypoint.R",
    "r-service/install.R",
    "r-service/Dockerfile",
    "backend/main.py",
    "backend/client.py",
    "backend/requirements.txt",
    "backend/Dockerfile",
    "frontend/Dockerfile",
    "docker-compose.yml",
]


@pytest.mark.parametrize("rel_path", REQUIRED_FILES)
def test_critical_files_exist(rel_path):
    assert (ROOT / rel_path).exists(), f"Missing required file: {rel_path}"


# ---------------------------------------------------------------------------
# 2. r-service port consistency
# ---------------------------------------------------------------------------

def test_dockerfile_port_matches_entrypoint():
    dockerfile = (ROOT / "r-service" / "Dockerfile").read_text()
    entrypoint = (ROOT / "r-service" / "entrypoint.R").read_text()

    expose_match = re.search(r"EXPOSE\s+(\d+)", dockerfile)
    assert expose_match, "No EXPOSE directive found in r-service/Dockerfile"
    expose_port = expose_match.group(1)

    # entrypoint.R: Sys.getenv("PORT", "8100")
    default_match = re.search(r'Sys\.getenv\("PORT",\s*"(\d+)"\)', entrypoint)
    assert default_match, 'No Sys.getenv("PORT", ...) found in r-service/entrypoint.R'
    default_port = default_match.group(1)

    assert expose_port == default_port, (
        f"Dockerfile EXPOSE {expose_port} != entrypoint.R default PORT {default_port}"
    )


# ---------------------------------------------------------------------------
# 3. plumber.R endpoint coverage
# ---------------------------------------------------------------------------

EXPECTED_ROUTES = [
    r"@post\s+/test-soft-powers",
    r"@post\s+/analyze",
    r"@post\s+/module-preservation",
    r"@post\s+/gene-selection",
    r"@get\s+/status/<job_id>",
    r"@get\s+/results/<job_id>/modules",
    r"@get\s+/results/<job_id>/plot",
    r"@get\s+/results/<job_id>/soft-power-plot",
    r"@get\s+/results/<job_id>/soft-powers",
    r"@get\s+/results/<job_id>/gene-selection",
    r"@get\s+/results/<job_id>/gene-selection-plot",
    r"@get\s+/results/<job_id>/preservation",
    r"@get\s+/results/<job_id>/preservation-plot",
    r"@get\s+/results/<job_id>/donor-counts",
]


@pytest.mark.parametrize("pattern", EXPECTED_ROUTES)
def test_plumber_endpoints_defined(pattern):
    plumber = (ROOT / "r-service" / "plumber.R").read_text()
    assert re.search(pattern, plumber), (
        f"plumber.R is missing an endpoint matching: {pattern}"
    )


def test_soft_power_table_is_extracted():
    """The numeric power table must be read out, not just plotted.

    Without GetPowerTable() the only soft-power output is a PNG and no caller can pick a
    power programmatically -- which is the bug this endpoint exists to fix.
    """
    plumber = (ROOT / "r-service" / "plumber.R").read_text()
    assert "GetPowerTable(" in plumber
    assert "soft_power_table.csv" in plumber


def test_gene_selection_uses_real_hdwgcna():
    """The sweep must call SetupForWGCNA/GetWGCNAGenes, never reimplement gene selection."""
    plumber = (ROOT / "r-service" / "plumber.R").read_text()
    assert "GetWGCNAGenes(" in plumber
    assert "gene_selection.csv" in plumber


def test_gene_select_not_hardcoded():
    """gene_select/fraction must come from the request, not baked into SetupForWGCNA."""
    plumber = (ROOT / "r-service" / "plumber.R").read_text()
    assert "gene_select = p$gene_select" in plumber
    assert "fraction    = p$fraction" in plumber


def test_soft_power_never_infinite():
    """ConstructNetwork(soft_power=NULL) resolves to Inf when no power clears 0.8.

    recommend_power() must return NA instead, and /analyze must refuse a non-finite power.
    """
    plumber = (ROOT / "r-service" / "plumber.R").read_text()
    assert "recommend_power" in plumber
    assert "NA_integer_" in plumber
    assert "is.finite(sp)" in plumber


# ---------------------------------------------------------------------------
# 4. docker-compose topology
# ---------------------------------------------------------------------------

def _compose() -> dict:
    return yaml.safe_load((ROOT / "docker-compose.yml").read_text())


def test_docker_compose_services_present():
    services = _compose()["services"]
    for name in ("r-service", "backend", "react-ui"):
        assert name in services, f"docker-compose.yml missing service: {name}"


def test_docker_compose_r_service_port():
    services = _compose()["services"]
    ports = services["r-service"].get("ports", [])
    assert any("8100" in str(p) for p in ports), (
        "r-service must expose port 8100"
    )


def test_docker_compose_backend_port():
    services = _compose()["services"]
    ports = services["backend"].get("ports", [])
    assert any("8200" in str(p) for p in ports), (
        "backend must expose port 8200"
    )


def test_docker_compose_shared_mount_on_both_services():
    # /shared is bind-mounted from the SHARED_DIR host directory, which lives OUTSIDE this
    # repo because the datasets are ~16.8 GB and must never enter git. The same folder is
    # visible at /shared in both containers on Windows and Linux, so users can drop .rds
    # files into it without docker cp. Both services must mount it, and the interpolation
    # must keep the ./local_data fallback so a bare `docker compose up` with no .env still
    # works. See .env.example and the README's "Where the data lives" section.
    services = _compose()["services"]
    for svc in ("r-service", "backend"):
        svc_volumes = [str(v) for v in services[svc].get("volumes", [])]
        shared = [v for v in svc_volumes if v.endswith(":/shared")]
        assert shared, f"service '{svc}' must bind-mount a host directory to /shared"
        assert any("${SHARED_DIR" in v for v in shared), (
            f"service '{svc}' must mount ${{SHARED_DIR}} at /shared, not a hardcoded path"
        )
        assert any(":-./local_data}" in v for v in shared), (
            f"service '{svc}' must keep the ./local_data default for SHARED_DIR"
        )


# ---------------------------------------------------------------------------
# 5. R_SERVICE_URL env var wiring
# ---------------------------------------------------------------------------

def test_client_reads_r_service_url_env_var():
    client_src = (ROOT / "backend" / "client.py").read_text()
    assert 'R_SERVICE_URL' in client_src, (
        "backend/client.py must read the R_SERVICE_URL environment variable"
    )


def test_docker_compose_sets_r_service_url():
    services = _compose()["services"]
    env = services["backend"].get("environment", {})
    # docker-compose env can be a list ("KEY=val") or a dict
    if isinstance(env, list):
        env_str = "\n".join(env)
    else:
        env_str = "\n".join(f"{k}={v}" for k, v in env.items())

    assert "R_SERVICE_URL" in env_str, (
        "docker-compose backend service must set R_SERVICE_URL"
    )
    assert "r-service" in env_str, (
        "R_SERVICE_URL must reference the 'r-service' hostname"
    )
    assert "8100" in env_str, (
        "R_SERVICE_URL must reference port 8100"
    )


# ---------------------------------------------------------------------------
# 6. Python syntax
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("src", list((ROOT / "backend").glob("*.py")))
def test_python_syntax(src):
    try:
        py_compile.compile(str(src), doraise=True)
    except py_compile.PyCompileError as e:
        pytest.fail(str(e))


# ---------------------------------------------------------------------------
# 7. Pydantic model — required-field enforcement
# ---------------------------------------------------------------------------

def test_pipeline_request_required_fields():
    from pydantic import ValidationError
    from main import PipelineRequest

    required = {
        "h5seurat_path": "/shared/data/foo.h5seurat",
        "out_dir": "/shared/results/foo",
        "cell_type_col": "cell_type",
        "group_by": ["sample", "cell_type"],
        "group_name": "Neuron",
        "wgcna_name": "test",
    }

    # Full valid request must not raise
    PipelineRequest(**required)

    # Dropping each required field must raise ValidationError
    for field in required:
        partial = {k: v for k, v in required.items() if k != field}
        with pytest.raises(ValidationError):
            PipelineRequest(**partial)


# ---------------------------------------------------------------------------
# 8. /analyze soft_power guard — unit level (no HTTP server)
# ---------------------------------------------------------------------------

def test_analyze_soft_power_guard():
    """
    The FastAPI route raises HTTPException(422) before forwarding to the R service
    when soft_power is absent. Verify the guard exists in the route handler source.
    """
    main_src = (ROOT / "backend" / "main.py").read_text()
    # Guard: `if body.soft_power is None`
    assert "soft_power is None" in main_src, (
        "backend/main.py /analyze route must guard against missing soft_power"
    )
    # Must raise an HTTPException, not just log or pass
    assert "HTTPException" in main_src, (
        "backend/main.py must import and raise HTTPException"
    )
