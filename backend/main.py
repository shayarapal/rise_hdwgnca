from fastapi import FastAPI, HTTPException, Response
from pydantic import BaseModel
from typing import Optional
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


@app.post("/test-soft-powers")
async def test_soft_powers(body: PipelineRequest):
    try:
        return await client.post_json("/test-soft-powers", body.model_dump(exclude_none=True))
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
