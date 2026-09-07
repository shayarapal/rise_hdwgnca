import { useState, useEffect, useCallback } from 'react'
import { listFiles, previewFile, fileDownloadUrl, triggerDownload } from '../api.js'

const fmtSize = (n) => {
  if (n == null) return ''
  if (n < 1024) return `${n} B`
  const u = ['KB', 'MB', 'GB', 'TB']
  let v = n / 1024, i = 0
  while (v >= 1024 && i < u.length - 1) { v /= 1024; i++ }
  return `${v.toFixed(1)} ${u[i]}`
}

const iconFor = (e) =>
  e.type === 'dir' ? '📁' :
  ['csv', 'tsv', 'tab', 'txt'].includes(e.ext) ? '📄' :
  ['png', 'jpg', 'jpeg', 'gif', 'svg'].includes(e.ext) ? '🖼️' :
  ['h5seurat', 'h5', 'rds', 'rdata'].includes(e.ext) ? '🧬' : '📃'

// mode: 'file' picks a file, 'dir' picks the current folder.
// browseOnly: read/preview/download only — no "use this" selection action.
export default function FileBrowser({ mode = 'file', title, onSelect, onClose, browseOnly = false }) {
  const [cwd, setCwd]         = useState(null)
  const [listing, setListing] = useState(null)
  const [selected, setSelected] = useState(null)  // selected file entry
  const [preview, setPreview] = useState(null)
  const [pathInput, setPathInput] = useState('')
  const [loading, setLoading] = useState(false)
  const [error, setError]     = useState(null)

  const load = useCallback(async (path) => {
    setLoading(true); setError(null)
    try {
      const data = await listFiles(path)
      setListing(data)
      setCwd(data.path)
      setPathInput(data.path)
      setSelected(null)
      setPreview(null)
    } catch (e) {
      setError(e.message)
    } finally {
      setLoading(false)
    }
  }, [])

  useEffect(() => { load(null) }, [load])

  const openEntry = useCallback(async (entry) => {
    if (entry.type === 'dir') { load(entry.path); return }
    setSelected(entry)
    setPreview({ loading: true })
    try {
      setPreview(await previewFile(entry.path))
    } catch (e) {
      setPreview({ kind: 'error', message: e.message })
    }
  }, [load])

  const chooseDisabled = mode === 'file' ? !selected : !cwd

  const handleChoose = () => {
    if (mode === 'file' && selected) onSelect(selected.path)
    else if (mode === 'dir' && cwd) onSelect(cwd)
  }

  return (
    <div className="fb-overlay" onClick={onClose}>
      <div className="fb-modal" onClick={(e) => e.stopPropagation()}>
        <div className="fb-header">
          <span>{title || (browseOnly ? 'Browse files' : mode === 'dir' ? 'Choose a folder' : 'Choose a file')}</span>
          <button className="fb-x" onClick={onClose} aria-label="Close">✕</button>
        </div>

        <div className="fb-toolbar">
          <button className="fb-btn" onClick={() => load(listing?.home)} title="Home">🏠</button>
          <button className="fb-btn" disabled={!listing?.parent} onClick={() => load(listing.parent)} title="Up">↑</button>
          <form
            className="fb-path"
            onSubmit={(e) => { e.preventDefault(); load(pathInput) }}
          >
            <input
              value={pathInput}
              onChange={(e) => setPathInput(e.target.value)}
              placeholder="/absolute/path"
              spellCheck={false}
            />
            <button className="fb-btn" type="submit">Go</button>
          </form>
        </div>

        <div className="fb-body">
          <div className="fb-list">
            {loading && <div className="fb-note">Loading…</div>}
            {error && <div className="fb-note fb-err">{error}</div>}
            {!loading && !error && listing?.entries.length === 0 && (
              <div className="fb-note">Empty folder</div>
            )}
            {!loading && !error && listing?.entries.map((e) => (
              <div
                key={e.path}
                className={`fb-row ${selected?.path === e.path ? 'sel' : ''}`}
                onClick={() => openEntry(e)}
                onDoubleClick={() => e.type === 'dir' && load(e.path)}
                title={e.path}
              >
                <span className="fb-icon">{iconFor(e)}</span>
                <span className="fb-name">{e.name}</span>
                <span className="fb-size">{fmtSize(e.size)}</span>
              </div>
            ))}
          </div>

          <div className="fb-preview">
            {!selected && <div className="fb-note">Select a file to preview it.</div>}
            {selected && <PreviewPane entry={selected} preview={preview} />}
          </div>
        </div>

        <div className="fb-footer">
          <span className="fb-cwd">{cwd}</span>
          <div className="fb-actions">
            {selected && (
              <button
                className="btn-sm"
                onClick={() => triggerDownload(fileDownloadUrl(selected.path), selected.name)}
              >
                Download
              </button>
            )}
            <button className="btn-sm" onClick={onClose}>{browseOnly ? 'Close' : 'Cancel'}</button>
            {!browseOnly && (
              <button className="btn-sm btn-sm-primary" disabled={chooseDisabled} onClick={handleChoose}>
                {mode === 'dir' ? 'Use this folder' : 'Use this file'}
              </button>
            )}
          </div>
        </div>
      </div>
    </div>
  )
}

function PreviewPane({ entry, preview }) {
  if (!preview || preview.loading) return <div className="fb-note">Reading {entry.name}…</div>
  if (preview.kind === 'error')  return <div className="fb-note fb-err">{preview.message}</div>

  if (preview.kind === 'image') {
    return <img className="fb-img" src={preview.data_uri} alt={entry.name} />
  }
  if (preview.kind === 'binary') {
    return <div className="fb-note">{preview.note || `Binary file (${fmtSize(preview.size)}) — no text preview.`} Use Download to export it.</div>
  }
  if (preview.kind === 'text') {
    return (
      <>
        {preview.truncated && <div className="fb-note fb-trunc">Showing the first part of a large file.</div>}
        <pre className="fb-text">{preview.content}</pre>
      </>
    )
  }
  if (preview.kind === 'table') {
    const d = preview.delimiter === '\t' ? 'TSV' : preview.delimiter === ',' ? 'CSV' : `delimiter "${preview.delimiter}"`
    return (
      <>
        <div className="fb-note">{d} · {preview.columns.length} columns · showing {preview.shown_rows} rows{preview.truncated ? ' (truncated)' : ''}</div>
        <div className="fb-table-wrap">
          <table className="fb-table">
            <thead>
              <tr>{preview.columns.map((c, i) => <th key={i}>{c}</th>)}</tr>
            </thead>
            <tbody>
              {preview.rows.map((r, ri) => (
                <tr key={ri}>{preview.columns.map((_, ci) => <td key={ci}>{r[ci] ?? ''}</td>)}</tr>
              ))}
            </tbody>
          </table>
        </div>
      </>
    )
  }
  return null
}
