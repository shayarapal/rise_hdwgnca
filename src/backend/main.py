from fastapi import FastAPI, HTTPException, Response
from pydantic import BaseModel, Field
from typing import Literal, Optional
import httpx

import client
from files import router as files_router

app = FastAPI(title="hdWGCNA Bridge Backend")
app.include_router(files_router)


class PipelineRequest(BaseModel):
    h5seurat_path: str
    out_dir: str
    cell_type_col: str
    group_by: list[str]
    group_name: str
    wgcna_name: str
    k: Optional[int] = 25
    max_shared: Optional[int] = 15
    network_type: Optional[str] = "signed"
    soft_power: Optional[int] = None  # required for /analyze, absent for /test-soft-powers

    # Condition split. Required for /module-preservation; on /test-soft-powers they are
    # optional and restrict the soft-power curve to the reference condition, so the power
    # is chosen from the network that will actually be built.
    condition_col: Optional[str] = None
    ref_group: Optional[str] = None
    query_group: Optional[str] = None
    n_permutations: Optional[int] = None   # r-service defaults to 250 (tutorial value)
    preservation_name: Optional[str] = None

    # Gene selection. Left None so the r-service owns the defaults ("fraction" / 0.05) and an
    # untouched request body stays byte-identical to what it sends today. Literal and the
    # numeric bounds give a 422 on bad input without hand-written validation; "custom" is a
    # legal SelectNetworkGenes value but the r-service rejects it, as there is no gene_list param.
    gene_select: Optional[Literal["variable", "fraction", "all", "custom"]] = None
    fraction: Optional[float] = Field(default=None, gt=0, le=1)
    fractions: Optional[list[float]] = None   # /gene-selection sweep override


@app.post("/test-soft-powers")
async def test_soft_powers(body: PipelineRequest):
    if body.condition_col is not None and body.ref_group is None:
        raise HTTPException(
            status_code=422, detail="ref_group is required when condition_col is given"
        )
    try:
        return await client.post_json("/test-soft-powers", body.model_dump(exclude_none=True))
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.post("/gene-selection")
async def gene_selection(body: PipelineRequest):
    if body.condition_col is not None and body.ref_group is None:
        raise HTTPException(
            status_code=422, detail="ref_group is required when condition_col is given"
        )
    try:
        return await client.post_json("/gene-selection", body.model_dump(exclude_none=True))
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.post("/analyze")
async def analyze(body: PipelineRequest):
    if body.soft_power is None:
        raise HTTPException(status_code=422, detail="soft_power is required for /analyze")
    try:
        return await client.post_json("/analyze", body.model_dump(exclude_none=True))
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.post("/module-preservation")
async def module_preservation(body: PipelineRequest):
    missing = [
        name
        for name in ("soft_power", "condition_col", "ref_group", "query_group")
        if getattr(body, name) is None
    ]
    if missing:
        raise HTTPException(
            status_code=422,
            detail=f"required for /module-preservation: {', '.join(missing)}",
        )
    if body.ref_group == body.query_group:
        raise HTTPException(status_code=422, detail="ref_group and query_group must differ")
    try:
        return await client.post_json("/module-preservation", body.model_dump(exclude_none=True))
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/status/{job_id}")
async def status(job_id: str):
    try:
        return await client.get_json(f"/status/{job_id}")
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/modules")
async def modules(job_id: str):
    try:
        return await client.get_json(f"/results/{job_id}/modules")
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/soft-powers")
async def soft_powers(job_id: str):
    try:
        return await client.get_json(f"/results/{job_id}/soft-powers")
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/gene-selection")
async def gene_selection_results(job_id: str):
    try:
        return await client.get_json(f"/results/{job_id}/gene-selection")
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/gene-selection-plot")
async def gene_selection_plot(job_id: str):
    try:
        data, content_type = await client.get_bytes(f"/results/{job_id}/gene-selection-plot")
        return Response(content=data, media_type=content_type)
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/preservation")
async def preservation(job_id: str):
    try:
        return await client.get_json(f"/results/{job_id}/preservation")
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/donor-counts")
async def donor_counts(job_id: str):
    try:
        return await client.get_json(f"/results/{job_id}/donor-counts")
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/plot")
async def plot(job_id: str):
    try:
        data, content_type = await client.get_bytes(f"/results/{job_id}/plot")
        return Response(content=data, media_type=content_type)
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/soft-power-plot")
async def soft_power_plot(job_id: str):
    try:
        data, content_type = await client.get_bytes(f"/results/{job_id}/soft-power-plot")
        return Response(content=data, media_type=content_type)
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)


@app.get("/results/{job_id}/preservation-plot")
async def preservation_plot(job_id: str):
    try:
        data, content_type = await client.get_bytes(f"/results/{job_id}/preservation-plot")
        return Response(content=data, media_type=content_type)
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=e.response.text)
