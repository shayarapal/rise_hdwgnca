"""
Phase 2 — Backend HTTP Smoke Tests
FastAPI layer tested via Starlette TestClient; the R service is mocked by
patching client.post_json / client.get_json / client.get_bytes directly.
No Docker or R process required.
"""
import unittest.mock as mock
from unittest.mock import AsyncMock, MagicMock

import httpx
import pytest
from starlette.testclient import TestClient

from main import app

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def tc():
    with TestClient(app) as c:
        yield c


VALID_TSP = {
    "h5seurat_path": "/shared/data/foo.h5seurat",
    "out_dir": "/shared/results/foo",
    "cell_type_col": "cell_type",
    "group_by": ["sample", "cell_type"],
    "group_name": "Neuron",
    "wgcna_name": "test",
}

VALID_ANALYZE = {**VALID_TSP, "soft_power": 6}

JOB_ID = "aaaaaaaa-0000-0000-0000-000000000001"


def _http_error(status: int) -> httpx.HTTPStatusError:
    mock_resp = MagicMock()
    mock_resp.status_code = status
    mock_resp.text = f"R service error {status}"
    return httpx.HTTPStatusError("error", request=MagicMock(), response=mock_resp)


# ---------------------------------------------------------------------------
# 1. Missing required fields → 422 (Pydantic, no R call)
# ---------------------------------------------------------------------------

def test_post_test_soft_powers_missing_fields(tc):
    # Empty body — all required fields absent
    r = tc.post("/test-soft-powers", json={})
    assert r.status_code == 422


def test_post_test_soft_powers_partial_fields(tc):
    # Only one required field present
    r = tc.post("/test-soft-powers", json={"h5seurat_path": "/shared/data/foo.h5seurat"})
    assert r.status_code == 422


# ---------------------------------------------------------------------------
# 2. /analyze requires soft_power
# ---------------------------------------------------------------------------

def test_post_analyze_missing_soft_power(tc):
    r = tc.post("/analyze", json=VALID_TSP)  # no soft_power
    assert r.status_code == 422


def test_post_analyze_with_soft_power_proxies(tc):
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/analyze", json=VALID_ANALYZE)
    assert r.status_code == 200
    assert r.json()["job_id"] == JOB_ID
    # Verify the body forwarded to R includes soft_power
    _, call_kwargs = m.call_args
    forwarded = m.call_args.args[1]  # second positional arg is the body dict
    assert forwarded["soft_power"] == 6


# ---------------------------------------------------------------------------
# 3. /test-soft-powers happy path proxies to R
# ---------------------------------------------------------------------------

def test_post_test_soft_powers_proxies(tc):
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/test-soft-powers", json=VALID_TSP)
    assert r.status_code == 200
    assert r.json()["job_id"] == JOB_ID
    m.assert_called_once()
    path_called = m.call_args.args[0]
    assert path_called == "/test-soft-powers"


# ---------------------------------------------------------------------------
# 4. /status proxies job_id through
# ---------------------------------------------------------------------------

def test_get_status_running(tc):
    with mock.patch("client.get_json", new=AsyncMock(return_value={"status": "running"})) as m:
        r = tc.get(f"/status/{JOB_ID}")
    assert r.status_code == 200
    assert r.json()["status"] == "running"
    assert JOB_ID in m.call_args.args[0]


def test_get_status_done(tc):
    with mock.patch("client.get_json", new=AsyncMock(return_value={"status": "done"})):
        r = tc.get(f"/status/{JOB_ID}")
    assert r.json()["status"] == "done"


def test_get_status_failed(tc):
    payload = {"status": "failed", "error": "LoadH5Seurat: file not found"}
    with mock.patch("client.get_json", new=AsyncMock(return_value=payload)):
        r = tc.get(f"/status/{JOB_ID}")
    assert r.json()["status"] == "failed"
    assert "error" in r.json()


# ---------------------------------------------------------------------------
# 5. /results/{job_id}/modules returns JSON
# ---------------------------------------------------------------------------

def test_get_modules_returns_json(tc):
    modules_payload = [
        {"gene": "APOE", "module": "M1", "kME": 0.87},
        {"gene": "CLU",  "module": "M1", "kME": 0.73},
    ]
    with mock.patch("client.get_json", new=AsyncMock(return_value=modules_payload)):
        r = tc.get(f"/results/{JOB_ID}/modules")
    assert r.status_code == 200
    assert isinstance(r.json(), list)
    assert r.json()[0]["gene"] == "APOE"


# ---------------------------------------------------------------------------
# 6. /results/{job_id}/plot and /soft-power-plot return image/png
# ---------------------------------------------------------------------------

_PNG_HEADER = b"\x89PNG\r\n\x1a\n" + b"\x00" * 8  # minimal PNG-like bytes


def test_get_network_plot_content_type(tc):
    with mock.patch("client.get_bytes", new=AsyncMock(return_value=(_PNG_HEADER, "image/png"))):
        r = tc.get(f"/results/{JOB_ID}/plot")
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"
    assert r.content == _PNG_HEADER


def test_get_soft_power_plot_content_type(tc):
    with mock.patch("client.get_bytes", new=AsyncMock(return_value=(_PNG_HEADER, "image/png"))):
        r = tc.get(f"/results/{JOB_ID}/soft-power-plot")
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"


# ---------------------------------------------------------------------------
# 7. R service errors propagate with correct status codes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("r_status", [400, 404, 500])
def test_r_service_error_propagates_post(tc, r_status):
    with mock.patch("client.post_json", new=AsyncMock(side_effect=_http_error(r_status))):
        r = tc.post("/test-soft-powers", json=VALID_TSP)
    assert r.status_code == r_status


@pytest.mark.parametrize("r_status", [404, 500])
def test_r_service_error_propagates_get_status(tc, r_status):
    with mock.patch("client.get_json", new=AsyncMock(side_effect=_http_error(r_status))):
        r = tc.get(f"/status/{JOB_ID}")
    assert r.status_code == r_status


@pytest.mark.parametrize("r_status", [404, 500])
def test_r_service_error_propagates_get_plot(tc, r_status):
    with mock.patch("client.get_bytes", new=AsyncMock(side_effect=_http_error(r_status))):
        r = tc.get(f"/results/{JOB_ID}/plot")
    assert r.status_code == r_status
