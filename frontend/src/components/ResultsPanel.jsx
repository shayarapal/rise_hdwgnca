import ModuleTable from './ModuleTable.jsx'
import { softPowerPlotUrl, networkPlotUrl, downloadTable, triggerDownload } from '../api.js'

function PlotBlock({ title, src, filename }) {
  return (
    <div className="result-block">
      <div className="result-block-header">
        <span>{title}</span>
        <button className="btn-sm" onClick={() => triggerDownload(src, filename)}>Export PNG</button>
      </div>
      <div className="result-block-body">
        <img className="result-img" src={src} alt={title} loading="lazy" />
      </div>
    </div>
  )
}

export default function ResultsPanel({ softPowerJobId, analyzeJobId, modules }) {
  const hasSoftPower = Boolean(softPowerJobId)
  const hasAnalyze   = Boolean(analyzeJobId && modules)

  if (!hasSoftPower) {
    return (
      <div className="empty-state">
        Run "Soft Power Test" to see results here.
      </div>
    )
  }

  return (
    <>
      <PlotBlock
        title="Soft Power Diagnostic"
        src={softPowerPlotUrl(softPowerJobId)}
        filename="soft-power-diagnostic.png"
      />

      {analyzeJobId && (
        <PlotBlock
          title="kME Plot (Module Hub Genes)"
          src={networkPlotUrl(analyzeJobId)}
          filename="kme-plot.png"
        />
      )}

      {modules && modules.length > 0 && (
        <div className="result-block">
          <div className="result-block-header">
            <span>Module Table — {modules.length} genes</span>
            <span className="result-header-actions">
              <button className="btn-sm" onClick={() => downloadTable(modules, 'modules.csv', ',')}>Export CSV</button>
              <button className="btn-sm" onClick={() => downloadTable(modules, 'modules.tsv', '\t')}>Export TSV</button>
            </span>
          </div>
          <div className="result-block-body" style={{ padding: 0 }}>
            <ModuleTable rows={modules} />
          </div>
        </div>
      )}
    </>
  )
}
