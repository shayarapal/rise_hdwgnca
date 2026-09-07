import { useState, useEffect, useCallback } from 'react'
import './styles.css'
import PipelineStepper from './components/PipelineStepper.jsx'
import ParamForm       from './components/ParamForm.jsx'
import JobStatus       from './components/JobStatus.jsx'
import ResultsPanel    from './components/ResultsPanel.jsx'
import FileBrowser     from './components/FileBrowser.jsx'
import {
  runSoftPowers, runAnalyze, runModulePreservation, runGeneSelection,
  getStatus, getModules, getPreservation, getSoftPowers, getGeneSelection,
} from './api.js'

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
  // Condition-comparison mode; null until the user enables it in ParamForm.
  condition_col:  null,
  ref_group:      null,
  query_group:    null,
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
  const [preservation,     setPreservation]     = useState(null)
  const [browseFiles,      setBrowseFiles]      = useState(false)
  const [softPowers,       setSoftPowers]       = useState(null)
  const [geneSelJob,       setGeneSelJob]       = useState(IDLE_JOB)
  const [geneSelection,    setGeneSelection]    = useState(null)

  // Comparison mode drives which endpoint the second phase calls.
  const compare = params.condition_col != null

  // ── Polling: soft power ────────────────────────────────────
  useEffect(() => {
    if (softPowerJob.status !== 'running') return
    const iv = setInterval(async () => {
      try {
        const { status, error } = await getStatus(softPowerJob.id)
        setSoftPowerJob((j) => ({ ...j, status, error: error ?? null }))
        if (status === 'done') {
          const sp = await getSoftPowers(softPowerJob.id)
          setSoftPowers(sp)
          // Prefill with the power ConstructNetwork would itself pick. Functional update so a
          // value the user already typed is never clobbered. Stays empty when nothing reaches
          // the scale-free threshold — SoftPowerTable shows the warning instead.
          if (sp.recommended_power != null) {
            setSoftPowerValue((v) => (v === '' ? String(sp.recommended_power) : v))
          }
        }
      } catch (e) {
        setSoftPowerJob((j) => ({ ...j, status: 'failed', error: e.message }))
      }
    }, 5000)
    return () => clearInterval(iv)
  }, [softPowerJob.status, softPowerJob.id])

  // ── Polling: gene-selection sweep ──────────────────────────
  // Independent of the soft-power job: sweeping does not invalidate a run, choosing a
  // different fraction does (see the invalidation effect below).
  useEffect(() => {
    if (geneSelJob.status !== 'running') return
    const iv = setInterval(async () => {
      try {
        const { status, error } = await getStatus(geneSelJob.id)
        setGeneSelJob((j) => ({ ...j, status, error: error ?? null }))
        if (status === 'done') {
          const rows = await getGeneSelection(geneSelJob.id)
          setGeneSelection(Array.isArray(rows) ? rows : [])
        }
      } catch (e) {
        setGeneSelJob((j) => ({ ...j, status: 'failed', error: e.message }))
      }
    }, 5000)
    return () => clearInterval(iv)
  }, [geneSelJob.status, geneSelJob.id])

  // ── Invalidation: gene set changed ─────────────────────────
  // gene_select/fraction change WHICH GENES the network is built from, which changes the
  // soft-power curve. A power chosen on the old gene set must not be carried into a run on a
  // new one, so the whole soft-power step resets.
  useEffect(() => {
    setSoftPowerJob(IDLE_JOB)
    setSoftPowers(null)
    setSoftPowerValue('')
    setAnalyzeJob(IDLE_JOB)
    setModules(null)
    setPreservation(null)
  }, [params.gene_select, params.fraction])

  // ── Polling: analyze ───────────────────────────────────────
  useEffect(() => {
    if (analyzeJob.status !== 'running') return
    const iv = setInterval(async () => {
      try {
        const { status, error } = await getStatus(analyzeJob.id)
        setAnalyzeJob((j) => ({ ...j, status, error: error ?? null }))
        if (status === 'done') {
          // Both /analyze and /module-preservation write modules.csv for the network
          // they built, so the module table is fetched the same way in either mode.
          const rows = await getModules(analyzeJob.id)
          setModules(Array.isArray(rows) ? rows : [])
          if (compare) {
            const pres = await getPreservation(analyzeJob.id)
            setPreservation(Array.isArray(pres) ? pres : [])
          }
        }
      } catch (e) {
        setAnalyzeJob((j) => ({ ...j, status: 'failed', error: e.message }))
      }
    }, 5000)
    return () => clearInterval(iv)
  }, [analyzeJob.status, analyzeJob.id, compare])

  // ── Actions ────────────────────────────────────────────────
  const handleRunSoftPowers = useCallback(async () => {
    setSoftPowerJob({ id: null, status: 'running', error: null })
    setAnalyzeJob(IDLE_JOB)
    setModules(null)
    setPreservation(null)
    setSoftPowerValue('')
    setSoftPowers(null)
    try {
      const { job_id } = await runSoftPowers(params)
      setSoftPowerJob({ id: job_id, status: 'running', error: null })
    } catch (e) {
      setSoftPowerJob({ id: null, status: 'failed', error: e.message })
    }
  }, [params])

  const handleRunGeneSelection = useCallback(async () => {
    setGeneSelJob({ id: null, status: 'running', error: null })
    setGeneSelection(null)
    try {
      const { job_id } = await runGeneSelection(params)
      setGeneSelJob({ id: job_id, status: 'running', error: null })
    } catch (e) {
      setGeneSelJob({ id: null, status: 'failed', error: e.message })
    }
  }, [params])

  // Adopting a fraction from the sweep goes through setParams, so the invalidation effect
  // above fires and the soft-power step resets — exactly as if it were typed by hand.
  const handleUseFraction = useCallback((fraction) => {
    setParams((p) => ({ ...p, gene_select: 'fraction', fraction }))
  }, [])

  const handleRunAnalyze = useCallback(async () => {
    setAnalyzeJob({ id: null, status: 'running', error: null })
    setModules(null)
    setPreservation(null)
    const body = { ...params, soft_power: Number(softPowerValue) }
    try {
      const { job_id } = compare
        ? await runModulePreservation(body)
        : await runAnalyze(body)
      setAnalyzeJob({ id: job_id, status: 'running', error: null })
    } catch (e) {
      setAnalyzeJob({ id: null, status: 'failed', error: e.message })
    }
  }, [params, softPowerValue, compare])

  const activeStep = deriveStep(softPowerJob.status, analyzeJob.status)

  // In comparison mode the soft-power curve is computed on the reference condition only,
  // so ref_group must be set before the first phase runs, not just the second.
  const conditionsReady = !compare || (
    params.condition_col.trim() &&
    (params.ref_group ?? '').trim() &&
    (params.query_group ?? '').trim() &&
    params.ref_group !== params.query_group
  )

  const canRunSoftPowers = (
    params.h5seurat_path.trim() &&
    params.out_dir.trim() &&
    params.cell_type_col.trim() &&
    params.group_name.trim() &&
    conditionsReady &&
    softPowerJob.status !== 'running'
  )

  const canRunAnalyze = (
    softPowerJob.status === 'done' &&
    softPowerValue.trim() !== '' &&
    !isNaN(Number(softPowerValue)) &&
    conditionsReady &&
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
          <ParamForm
            params={params}
            onChange={setParams}
            onRunGeneSelection={handleRunGeneSelection}
            geneSelectionRunning={geneSelJob.status === 'running'}
          />

          <JobStatus
            status={geneSelJob.status}
            jobId={geneSelJob.id}
            error={geneSelJob.error}
            label="Gene selection"
          />

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
                <span className="hint">
                  {softPowers?.recommended_power != null
                    ? `Prefilled with ${softPowers.recommended_power} — the lowest power reaching R² ≥ ${softPowers.sft_threshold}. See the table on the right.`
                    : 'Pick the lowest power where scale-free fit ≥ 0.80'}
                </span>
              </div>

              <button
                className="btn btn-warn"
                onClick={handleRunAnalyze}
                disabled={!canRunAnalyze}
              >
                {compare ? 'Run Module Preservation' : 'Run Full Pipeline'}
              </button>

              <JobStatus
                status={analyzeJob.status}
                jobId={analyzeJob.id}
                error={analyzeJob.error}
                label={compare ? 'Module preservation' : 'Full pipeline'}
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
            preservation={preservation}
            refGroup={params.ref_group}
            queryGroup={params.query_group}
            softPowers={softPowers}
            geneSelectionJobId={geneSelJob.status === 'done' ? geneSelJob.id : null}
            geneSelection={geneSelection}
            onUseFraction={handleUseFraction}
            outDir={params.out_dir}
          />
        </section>
      </div>
    </>
  )
}
