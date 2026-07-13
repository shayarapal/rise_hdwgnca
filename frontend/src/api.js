// All paths are relative — Vite dev proxy forwards them to the backend on :8200

export async function runSoftPowers(params) {
  const r = await fetch('/test-soft-powers', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!r.ok) throw new Error(await r.text());
  return r.json(); // { job_id }
}

export async function runAnalyze(params) {
  const r = await fetch('/analyze', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!r.ok) throw new Error(await r.text());
  return r.json(); // { job_id }
}

export async function getStatus(jobId) {
  const r = await fetch(`/status/${jobId}`);
  if (!r.ok) throw new Error(await r.text());
  return r.json(); // { status: 'running'|'done'|'failed', error? }
}

export async function getModules(jobId) {
  const r = await fetch(`/results/${jobId}/modules`);
  if (!r.ok) throw new Error(await r.text());
  return r.json(); // array of module row objects
}

// Image URLs — used directly as <img src={...}> to avoid loading into JS memory
export const softPowerPlotUrl = (jobId) => `/results/${jobId}/soft-power-plot`;
export const networkPlotUrl   = (jobId) => `/results/${jobId}/plot`;

// ── Local filesystem browse / preview / export ────────────────────────────────

export async function listFiles(path) {
  const q = path ? `?path=${encodeURIComponent(path)}` : '';
  const r = await fetch(`/files/list${q}`);
  if (!r.ok) throw new Error(await r.text());
  return r.json(); // { path, parent, home, entries: [{ name, path, type, size, ext, modified }] }
}

export async function previewFile(path, rows = 200) {
  const r = await fetch(`/files/preview?path=${encodeURIComponent(path)}&rows=${rows}`);
  if (!r.ok) throw new Error(await r.text());
  return r.json(); // { kind: 'table'|'text'|'image'|'binary', ... }
}

// Direct URL for downloading (exporting) any file on the machine to the browser.
export const fileDownloadUrl = (path) => `/files/download?path=${encodeURIComponent(path)}`;

// Client-side export: turn in-memory rows into a downloaded csv/tsv file.
export function downloadTable(rows, filename, delimiter = ',') {
  if (!rows || rows.length === 0) return;
  const cols = Object.keys(rows[0]);
  const esc = (v) => {
    const s = v == null ? '' : String(v);
    return /["\n\r]|[,\t]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  const lines = [cols.join(delimiter)];
  for (const row of rows) lines.push(cols.map((c) => esc(row[c])).join(delimiter));
  const blob = new Blob([lines.join('\n')], { type: 'text/plain;charset=utf-8' });
  triggerDownload(URL.createObjectURL(blob), filename, true);
}

// Trigger a browser download for a URL (blob or server route).
export function triggerDownload(url, filename, revoke = false) {
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  if (revoke) setTimeout(() => URL.revokeObjectURL(url), 1000);
}
