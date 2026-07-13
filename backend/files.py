# Local filesystem browse / preview / download.
#
# Running locally, the backend process shares the user's filesystem, so these
# endpoints let the web UI read files and folders from the machine (to pick an
# h5Seurat input, an output dir, or inspect a csv/tsv) and export any file back
# to the browser. Under Docker the backend only sees the /shared volume.
import base64
import csv
import io
import mimetypes
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse

router = APIRouter(prefix="/files", tags=["files"])

HOME = str(Path.home())

# Extensions we force a tabular delimiter for; everything else is sniffed.
_DELIMS = {"csv": ",", "tsv": "\t", "tab": "\t"}
# Cap how much of a file we read into memory for a preview.
_TEXT_BYTES = 1_000_000
_IMG_BYTES = 5 * 1024 * 1024


def _resolve(path: str) -> Path:
    try:
        return Path(path).expanduser().resolve()
    except (OSError, RuntimeError, ValueError):
        raise HTTPException(status_code=400, detail=f"Invalid path: {path}")


@router.get("/list")
def list_dir(path: Optional[str] = Query(default=None)):
    """List folders and files at `path` (defaults to the user's home dir)."""
    p = _resolve(path) if path else Path(HOME)
    if not p.exists():
        raise HTTPException(status_code=404, detail=f"Not found: {p}")
    if not p.is_dir():
        raise HTTPException(status_code=400, detail=f"Not a directory: {p}")

    entries = []
    try:
        children = list(p.iterdir())
    except PermissionError:
        raise HTTPException(status_code=403, detail=f"Permission denied: {p}")

    # Folders first, then files, each alphabetical (case-insensitive).
    for child in sorted(children, key=lambda c: (c.is_file(), c.name.lower())):
        try:
            st = child.stat()
            is_dir = child.is_dir()
            entries.append({
                "name": child.name,
                "path": str(child),
                "type": "dir" if is_dir else "file",
                "size": None if is_dir else st.st_size,
                "ext": child.suffix.lower().lstrip("."),
                "modified": st.st_mtime,
            })
        except (PermissionError, OSError):
            continue  # skip unreadable entries rather than failing the listing

    return {
        "path": str(p),
        "parent": str(p.parent) if p.parent != p else None,
        "home": HOME,
        "entries": entries,
    }


@router.get("/preview")
def preview(path: str, rows: int = Query(default=200, ge=1, le=5000)):
    """Read any file for display: tabular (csv/tsv/sniffed), text, image, or binary."""
    p = _resolve(path)
    if not p.exists() or not p.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {p}")

    size = p.stat().st_size
    ext = p.suffix.lower().lstrip(".")
    mime, _ = mimetypes.guess_type(str(p))

    if mime and mime.startswith("image/"):
        if size > _IMG_BYTES:
            # Too large to inline as a data URI; a truncated image won't render.
            return {"kind": "binary", "size": size, "ext": ext, "mime": mime,
                    "note": "Image too large to preview inline — use Download."}
        with p.open("rb") as f:  # bounded read; never load the whole file
            data = f.read(_IMG_BYTES)
        b64 = base64.b64encode(data).decode("ascii")
        return {"kind": "image", "mime": mime, "size": size,
                "data_uri": f"data:{mime};base64,{b64}"}

    try:
        with p.open("r", encoding="utf-8", errors="replace") as f:
            head = f.read(_TEXT_BYTES)
    except (OSError, UnicodeError):
        return {"kind": "binary", "size": size, "ext": ext, "mime": mime}

    if "\x00" in head:  # NUL byte → treat as binary
        return {"kind": "binary", "size": size, "ext": ext, "mime": mime}

    delim = _DELIMS.get(ext)
    if delim is None:
        try:
            sample = "\n".join(head.splitlines()[:20])
            delim = csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter
        except (csv.Error, IndexError):
            delim = None

    if delim:
        parsed = list(csv.reader(io.StringIO(head), delimiter=delim))
        header = parsed[0] if parsed else []
        body = parsed[1:rows + 1]
        return {
            "kind": "table",
            "delimiter": delim,
            "columns": header,
            "rows": body,
            "shown_rows": len(body),
            "size": size,
            "truncated": size > len(head.encode("utf-8", "replace")),
        }

    return {
        "kind": "text",
        "content": head,
        "size": size,
        "truncated": size > len(head.encode("utf-8", "replace")),
    }


@router.get("/download")
def download(path: str):
    """Stream a file to the browser so the user can export it to their computer."""
    p = _resolve(path)
    if not p.exists() or not p.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {p}")
    return FileResponse(str(p), filename=p.name, media_type="application/octet-stream")
