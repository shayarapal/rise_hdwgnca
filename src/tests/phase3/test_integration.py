"""
Phase 3 — Docker Compose Integration Tests
Requires `make docker-up` to be running before executing this suite.
All pipeline jobs submitted here use a fake h5seurat path and are expected
to fail inside R (file not found) — that is the correct observable behaviour.
"""
import pytest

from .conftest import R_SERVICE, BACKEND, poll_until_terminal

# ---------------------------------------------------------------------------
# Shared test payload — fake file path intentionally causes R-side failure
# ---------------------------------------------------------------------------

FAKE_PAYLOAD = {
    "h5seurat_path": "/shared/data/__nonexistent__.h5seurat",
    "out_dir":       "/shared/results/__phase3_test__",
    "cell_type_col": "cell_type",
    "group_by":      ["sample", "cell_type"],
    "group_name":    "Neuron",
    "wgcna_name":    "phase3_test",
}


# ---------------------------------------------------------------------------
# 1. Reachability
# ---------------------------------------------------------------------------

def test_r_service_reachable(r_client):
    """Any HTTP response (even 404) proves the R service is up."""
    try:
        r_client.get("/")
    except Exception as exc:
        pytest.fail(f"r-service not reachable: {exc}")


def test_backend_reachable(backend):
    try:
        backend.get("/")
    except Exception as exc:
        pytest.fail(f"backend not reachable: {exc}")


# ---------------------------------------------------------------------------
# 2. Backend field validation (no R involvement)
# ---------------------------------------------------------------------------

def test_backend_missing_fields_422(backend):
    r = backend.post("/test-soft-powers", json={})
    assert r.status_code == 422, r.text


def test_backend_analyze_missing_soft_power_422(backend):
    r = backend.post("/analyze", json=FAKE_PAYLOAD)
    assert r.status_code == 422, r.text


# ---------------------------------------------------------------------------
# 3. R service field validation (direct hit, bypassing backend)
# ---------------------------------------------------------------------------

def test_r_service_missing_fields_400(r_client):
    r = r_client.post("/test-soft-powers", json={})
    assert r.status_code == 400, r.text
    body = r.json()
    assert "error" in body
    assert "missing" in body["error"].lower()


# ---------------------------------------------------------------------------
# 4. Unknown job → 404
# ---------------------------------------------------------------------------

def test_unknown_job_status_404(backend):
    r = backend.get("/status/does-not-exist-00000000")
    assert r.status_code == 404, r.text


def test_unknown_job_modules_404(backend):
    r = backend.get("/results/does-not-exist-00000000/modules")
    assert r.status_code == 404, r.text


def test_unknown_job_plot_404(backend):
    r = backend.get("/results/does-not-exist-00000000/plot")
    assert r.status_code == 404, r.text


def test_unknown_job_soft_power_plot_404(backend):
    r = backend.get("/results/does-not-exist-00000000/soft-power-plot")
    assert r.status_code == 404, r.text


# ---------------------------------------------------------------------------
# 5. Full job lifecycle with fake data (expected to fail inside R)
# ---------------------------------------------------------------------------

def test_fake_path_job_submitted_and_fails(backend):
    """
    Submit a job with a non-existent h5seurat path.
    The job_id must be returned immediately, and when polled the job
    must eventually reach 'failed' with a non-empty error message.
    """
    payload = {**FAKE_PAYLOAD, "soft_power": 6}
    r = backend.post("/analyze", json=payload)
    assert r.status_code == 200, r.text

    body = r.json()
    assert "job_id" in body, f"Expected job_id in response, got: {body}"
    job_id = body["job_id"]
    assert job_id  # non-empty UUID

    terminal = poll_until_terminal(backend, job_id, timeout=90)

    assert terminal["status"] == "failed", (
        f"Expected status=failed for a fake file path, got: {terminal}"
    )
    assert terminal.get("error"), (
        "Expected a non-empty error message for a failed job"
    )


def test_fake_path_test_soft_powers_fails(backend):
    """Same lifecycle check for /test-soft-powers."""
    r = backend.post("/test-soft-powers", json=FAKE_PAYLOAD)
    assert r.status_code == 200, r.text

    job_id = r.json()["job_id"]
    terminal = poll_until_terminal(backend, job_id, timeout=90)

    assert terminal["status"] == "failed"
    assert terminal.get("error")
