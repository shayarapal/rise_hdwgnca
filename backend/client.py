import os
import httpx

_BASE = os.getenv("R_SERVICE_URL", "http://localhost:8100").rstrip("/")
_TIMEOUT = httpx.Timeout(10.0, read=None)  # long read for large result payloads


async def post_json(path: str, body: dict) -> dict:
    async with httpx.AsyncClient(base_url=_BASE, timeout=_TIMEOUT) as c:
        r = await c.post(path, json=body)
        r.raise_for_status()
        return r.json()


async def get_json(path: str) -> dict:
    async with httpx.AsyncClient(base_url=_BASE, timeout=_TIMEOUT) as c:
        r = await c.get(path)
        r.raise_for_status()
        return r.json()


async def get_bytes(path: str) -> tuple[bytes, str]:
    """Returns (content_bytes, content_type)."""
    async with httpx.AsyncClient(base_url=_BASE, timeout=_TIMEOUT) as c:
        r = await c.get(path)
        r.raise_for_status()
        return r.content, r.headers.get("content-type", "application/octet-stream")
