"""
Phase 3 session fixtures.
Tests are skipped automatically when the Docker stack isn't running.
Start the stack with:  make docker-up
"""
import time

import httpx
import pytest

R_SERVICE = "http://localhost:8100"
BACKEND   = "http://localhost:8200"


def _reachable(url: str) -> bool:
    try:
        httpx.get(url, timeout=3.0)
        return True
    except Exception:
        return False


@pytest.fixture(scope="session")
def services_up():
    r_ok = _reachable(R_SERVICE)
    b_ok = _reachable(BACKEND)
    if not r_ok or not b_ok:
        status = (
            f"r-service={'UP' if r_ok else 'DOWN'}, "
            f"backend={'UP' if b_ok else 'DOWN'}"
        )
        pytest.skip(f"Docker stack not running ({status}). Run 'make docker-up' first.")


@pytest.fixture(scope="session")
def backend(services_up):
    return httpx.Client(base_url=BACKEND, timeout=httpx.Timeout(10.0, read=None))


@pytest.fixture(scope="session")
def r_client(services_up):
    return httpx.Client(base_url=R_SERVICE, timeout=httpx.Timeout(10.0, read=None))


def poll_until_terminal(backend_client, job_id: str, timeout: int = 60) -> dict:
    """Poll /status/{job_id} until status is 'done' or 'failed', then return the status dict."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        r = backend_client.get(f"/status/{job_id}")
        r.raise_for_status()
        body = r.json()
        if body.get("status") in ("done", "failed"):
            return body
        time.sleep(2)
    raise TimeoutError(f"job {job_id} did not reach a terminal state within {timeout}s")
