import { useState, useEffect, useCallback } from 'react'
import './styles.css'
import PipelineStepper from './components/PipelineStepper.jsx'
import ParamForm       from './components/ParamForm.jsx'
import JobStatus       from './components/JobStatus.jsx'
import ResultsPanel    from './components/ResultsPanel.jsx'
import FileBrowser     from './components/FileBrowser.jsx'
import { runSoftPowers, runAnalyze, getStatus, getModules } from './api.js'

const DEFAULT_PARAMS = {
  h5seurat_path: '',
  out_dir:        '',
  cell_type_col:  'cell_type',
  group_by:       ['cell_type', 'sample_id'],
  group_name:     '',
  wgcna_name:     'tutorial',
  k:              25,
  max_shared:     15,
  network_type:   'signed',
}

const IDLE_JOB = { id: null, status: 'idle', error: null }

// Derives the active pipeline step (0–5) from job states
function deriveStep(spStatus, anStatus) {
  if (anStatus === 'done')    return 5
  if (anStatus === 'running') return 4
  if (spStatus === 'done')    return 3
  if (spStatus === 'running') return 2
  return 0
}

export default function App() {
  const [params,           setParams]          = useState(DEFAULT_PARAMS)
  const [softPowerJob,     setSoftPowerJob]     = useState(IDLE_JOB)
  const [analyzeJob,       setAnalyzeJob]       = useState(IDLE_JOB)
  const [softPowerValue,   setSoftPowerValue]   = useState('')
  const [modules,          setModules]          = useState(null)
  const [browseFiles,      setBrowseFiles]      = useState(false)

  // ── Polling: soft power ────────────────────────────────────
  useEffect(() => {
    if (softPowerJob.status !== 'running') return
    const iv = setInterval(async () => {
      try {
        const { status, error } = await getStatus(softPowerJob.id)
        setSoftPowerJob((j) => ({ ...j, status, error: error ?? null }))
      } catch (e) {
        setSoftPowerJob((j) => ({ ...j, status: 'failed', error: e.message }))
      }
    }, 5000)
    return () => clearInterval(iv)
  }, [softPowerJob.status, softPowerJob.id])

  // ── Polling: analyze ───────────────────────────────────────
  useEffect(() => {
    if (analyzeJob.status !== 'running') return
    const iv = setInterval(async () => {
      try {
        const { status, error } = await getStatus(analyzeJob.id)
        setAnalyzeJob((j) => ({ ...j, status, error: error ?? null }))
        if (status === 'done') {
          const rows = await getModules(analyzeJob.id)
          setModules(Array.isArray(rows) ? rows : [])
        }
      } catch (e) {
        setAnalyzeJob((j) => ({ ...j, status: 'failed', error: e.message }))
      }
    }, 5000)
    return () => clearInterval(iv)
  }, [analyzeJob.status, analyzeJob.id])

  // ── Actions ────────────────────────────────────────────────
  const handleRunSoftPowers = useCallback(async () => {
    setSoftPowerJob({ id: null, status: 'running', error: null })
    setAnalyzeJob(IDLE_JOB)
    setModules(null)
    setSoftPowerValue('')
    try {
      const { job_id } = await runSoftPowers(params)
      setSoftPowerJob({ id: job_id, status: 'running', error: null })
    } catch (e) {
      setSoftPowerJob({ id: null, status: 'failed', error: e.message })
    }
  }, [params])

  const handleRunAnalyze = useCallback(async () => {
    setAnalyzeJob({ id: null, status: 'running', error: null })
    setModules(null)
    try {
      const { job_id } = await runAnalyze({
        ...params,
        soft_power: Number(softPowerValue),
      })
      setAnalyzeJob({ id: job_id, status: 'running', error: null })
    } catch (e) {
      setAnalyzeJob({ id: null, status: 'failed', error: e.message })
    }
  }, [params, softPowerValue])

  const activeStep = deriveStep(softPowerJob.status, analyzeJob.status)

  const canRunSoftPowers = (
    params.h5seurat_path.trim() &&
    params.out_dir.trim() &&
    params.cell_type_col.trim() &&
    params.group_name.trim() &&
    softPowerJob.status !== 'running'
  )

  const canRunAnalyze = (
    softPowerJob.status === 'done' &&
    softPowerValue.trim() !== '' &&
    !isNaN(Number(softPowerValue)) &&
    analyzeJob.status !== 'running'
  )

  return (
    <>
      <header className="app-header">
        <h1>hdWGCNA Bridge</h1>
        <span className="subtitle">Single-cell co-expression network analysis</span>
        <button className="btn-header" onClick={() => setBrowseFiles(true)}>📂 Browse files</button>
      </header>

      {browseFiles && (
        <FileBrowser browseOnly title="Browse files on this computer" onClose={() => setBrowseFiles(false)} />
      )}

      <PipelineStepper activeStep={activeStep} />

      <div className="workspace">
        {/* ── Left: Parameters ── */}
        <aside className="params-panel">
          <ParamForm params={params} onChange={setParams} />

          <button
            className="btn btn-primary"
            onClick={handleRunSoftPowers}
            disabled={!canRunSoftPowers}
          >
            Run Soft Power Test
          </button>

          <JobStatus
            status={softPowerJob.status}
            jobId={softPowerJob.id}
            error={softPowerJob.error}
            label="Soft power"
          />

          {/* ── Break point: soft power selection ── */}
          {softPowerJob.status === 'done' && (
            <>
              <div className="break-divider">Inspect the plot → choose soft power</div>

              <div className="form-group">
                <label htmlFor="soft_power">Soft power threshold</label>
                <input
                  id="soft_power"
                  type="number"
                  min={1}
                  max={30}
                  value={softPowerValue}
                  onChange={(e) => setSoftPowerValue(e.target.value)}
                  placeholder="e.g. 9"
                  autoFocus
                />
                <span className="hint">Pick the lowest power where scale-free fit ≥ 0.80</span>
              </div>

              <button
                className="btn btn-warn"
                onClick={handleRunAnalyze}
                disabled={!canRunAnalyze}
              >
                Run Full Pipeline
              </button>

              <JobStatus
                status={analyzeJob.status}
                jobId={analyzeJob.id}
                error={analyzeJob.error}
                label="Full pipeline"
              />
            </>
          )}
        </aside>

        {/* ── Right: Results ── */}
        <section className="results-panel">
          <ResultsPanel
            softPowerJobId={softPowerJob.status === 'done' ? softPowerJob.id : null}
            analyzeJobId={analyzeJob.status === 'done' ? analyzeJob.id : null}
            modules={modules}
          />
        </section>
      </div>
    </>
  )
}
