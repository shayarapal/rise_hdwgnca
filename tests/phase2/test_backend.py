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

VALID_PRESERVATION = {
    **VALID_ANALYZE,
    "condition_col": "condition",
    "ref_group": "control",
    "query_group": "PD",
}

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
# 3b. /module-preservation required fields
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dropped", ["soft_power", "condition_col", "ref_group", "query_group"])
def test_post_module_preservation_missing_field(tc, dropped):
    body = {k: v for k, v in VALID_PRESERVATION.items() if k != dropped}
    r = tc.post("/module-preservation", json=body)
    assert r.status_code == 422
    assert dropped in r.json()["detail"]


def test_post_module_preservation_same_ref_and_query(tc):
    body = {**VALID_PRESERVATION, "query_group": "control"}
    r = tc.post("/module-preservation", json=body)
    assert r.status_code == 422


def test_post_module_preservation_proxies(tc):
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/module-preservation", json=VALID_PRESERVATION)
    assert r.status_code == 200
    assert r.json()["job_id"] == JOB_ID
    assert m.call_args.args[0] == "/module-preservation"
    forwarded = m.call_args.args[1]
    assert forwarded["ref_group"] == "control"
    assert forwarded["query_group"] == "PD"
    # n_permutations is left unset so the r-service applies the tutorial default (250)
    assert "n_permutations" not in forwarded


# ---------------------------------------------------------------------------
# 3c. Gene selection: params + sweep
# ---------------------------------------------------------------------------

def test_gene_select_and_fraction_omitted_when_unset(tc):
    """An untouched request must not carry gene_select/fraction, so the r-service applies
    its own defaults ("fraction" / 0.05) and the wire body stays what it is today."""
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/test-soft-powers", json=VALID_TSP)
    assert r.status_code == 200
    forwarded = m.call_args.args[1]
    assert "gene_select" not in forwarded
    assert "fraction" not in forwarded


def test_gene_select_and_fraction_forwarded(tc):
    body = {**VALID_TSP, "gene_select": "fraction", "fraction": 0.1}
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/test-soft-powers", json=body)
    assert r.status_code == 200
    forwarded = m.call_args.args[1]
    assert forwarded["gene_select"] == "fraction"
    assert forwarded["fraction"] == 0.1


def test_bad_gene_select_rejected(tc):
    r = tc.post("/test-soft-powers", json={**VALID_TSP, "gene_select": "bogus"})
    assert r.status_code == 422


@pytest.mark.parametrize("bad", [0, 1.5, -0.1])
def test_fraction_out_of_range_rejected(tc, bad):
    r = tc.post("/test-soft-powers", json={**VALID_TSP, "fraction": bad})
    assert r.status_code == 422


def test_post_gene_selection_proxies(tc):
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/gene-selection", json=VALID_TSP)
    assert r.status_code == 200
    assert r.json()["job_id"] == JOB_ID
    assert m.call_args.args[0] == "/gene-selection"


def test_post_gene_selection_forwards_fractions(tc):
    body = {**VALID_TSP, "fractions": [0.05, 0.2]}
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/gene-selection", json=body)
    assert r.status_code == 200
    assert m.call_args.args[1]["fractions"] == [0.05, 0.2]


def test_get_soft_powers_returns_table_and_recommendations(tc):
    payload = {
        "table": [{"Power": 1, "SFT.R.sq": 0.02, "slope": 0.9, "mean.k.": 900.0}],
        "recommended_power": 6,
        "smallest_power": 2,
        "differ": True,
        "sft_threshold": 0.8,
        "min_power": 3,
        "max_sft_r_sq": 0.91,
        "warning": None,
    }
    with mock.patch("client.get_json", new=AsyncMock(return_value=payload)):
        r = tc.get(f"/results/{JOB_ID}/soft-powers")
    assert r.status_code == 200
    body = r.json()
    assert body["recommended_power"] == 6
    assert body["smallest_power"] == 2
    assert body["table"][0]["SFT.R.sq"] == 0.02


def test_get_soft_powers_null_when_no_power_qualifies(tc):
    """The Inf bug: when nothing reaches 0.8 the recommendation must be null, not a number."""
    payload = {
        "table": [{"Power": 1, "SFT.R.sq": 0.1}],
        "recommended_power": None,
        "smallest_power": None,
        "differ": False,
        "sft_threshold": 0.8,
        "min_power": 3,
        "max_sft_r_sq": 0.62,
        "warning": "no power in the tested grid reaches SFT.R.sq >= 0.8",
    }
    with mock.patch("client.get_json", new=AsyncMock(return_value=payload)):
        r = tc.get(f"/results/{JOB_ID}/soft-powers")
    assert r.status_code == 200
    assert r.json()["recommended_power"] is None
    assert r.json()["warning"]


def test_get_gene_selection_returns_rows(tc):
    rows = [{"fraction": 0.05, "n_genes": 7465, "n_all_genes": 21455,
             "pct_of_all_genes": 34.79, "warning": None, "error": None}]
    with mock.patch("client.get_json", new=AsyncMock(return_value=rows)):
        r = tc.get(f"/results/{JOB_ID}/gene-selection")
    assert r.status_code == 200
    assert r.json()[0]["n_genes"] == 7465


def test_get_gene_selection_plot_is_png(tc):
    with mock.patch("client.get_bytes", new=AsyncMock(return_value=(b"\x89PNG_x", "image/png"))):
        r = tc.get(f"/results/{JOB_ID}/gene-selection-plot")
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"


# /test-soft-powers takes condition_col optionally, but not without ref_group.
def test_post_test_soft_powers_condition_col_without_ref_group(tc):
    body = {**VALID_TSP, "condition_col": "condition"}
    r = tc.post("/test-soft-powers", json=body)
    assert r.status_code == 422


def test_post_test_soft_powers_with_condition_forwards_ref_group(tc):
    body = {**VALID_TSP, "condition_col": "condition", "ref_group": "control"}
    with mock.patch("client.post_json", new=AsyncMock(return_value={"job_id": JOB_ID})) as m:
        r = tc.post("/test-soft-powers", json=body)
    assert r.status_code == 200
    forwarded = m.call_args.args[1]
    assert forwarded["condition_col"] == "condition"
    assert forwarded["ref_group"] == "control"


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


def test_get_preservation_plot_content_type(tc):
    with mock.patch("client.get_bytes", new=AsyncMock(return_value=(_PNG_HEADER, "image/png"))):
        r = tc.get(f"/results/{JOB_ID}/preservation-plot")
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"


# ---------------------------------------------------------------------------
# 6b. /results/{job_id}/preservation and /donor-counts return JSON
# ---------------------------------------------------------------------------

def test_get_preservation_returns_json(tc):
    payload = [
        {"module": "turquoise", "Zsummary.pres": 1.4, "moduleSize.Z": 812},
        {"module": "blue", "Zsummary.pres": 14.9, "moduleSize.Z": 402},
    ]
    with mock.patch("client.get_json", new=AsyncMock(return_value=payload)):
        r = tc.get(f"/results/{JOB_ID}/preservation")
    assert r.status_code == 200
    assert r.json()[0]["module"] == "turquoise"


def test_get_donor_counts_returns_json(tc):
    payload = [{"sample_id": "s_0096", "condition": "PD", "n_cells": 412,
                "clears_min_cells": True}]
    with mock.patch("client.get_json", new=AsyncMock(return_value=payload)):
        r = tc.get(f"/results/{JOB_ID}/donor-counts")
    assert r.status_code == 200
    assert r.json()[0]["clears_min_cells"] is True


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
