const LABELS = {
  idle:    '',
  running: 'Running…',
  done:    'Done',
  failed:  'Failed',
}

const ICONS = {
  running: <span className="spinner" aria-hidden="true" />,
  done:    '✓',
  failed:  '✕',
}

export default function JobStatus({ status, jobId, error, label }) {
  if (status === 'idle') return null

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
      <div className={`job-status ${status}`} role="status" aria-live="polite">
        {ICONS[status]}
        <span>{label ? `${label}: ` : ''}{LABELS[status]}</span>
      </div>
      {jobId && (
        <p className="job-id-text">job: {jobId}</p>
      )}
      {status === 'failed' && error && (
        <p className="error-text">{error}</p>
      )}
    </div>
  )
}
