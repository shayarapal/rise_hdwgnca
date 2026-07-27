import ModuleTable from './ModuleTable.jsx'
import PreservationTable from './PreservationTable.jsx'
import SoftPowerTable from './SoftPowerTable.jsx'
import GeneSelectionPanel from './GeneSelectionPanel.jsx'
import {
  softPowerPlotUrl, networkPlotUrl, preservationPlotUrl,
  downloadTable, triggerDownload, fileDownloadUrl,
} from '../api.js'

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

export default function ResultsPanel({
  softPowerJobId, analyzeJobId, modules, preservation, refGroup, queryGroup,
  softPowers, geneSelectionJobId, geneSelection, onUseFraction, outDir,
}) {
  const hasSoftPower = Boolean(softPowerJobId)
  const hasPreservation = Boolean(preservation && preservation.length > 0)
  const hasGeneSelection = Boolean(geneSelectionJobId && geneSelection?.length)

  if (!hasSoftPower && !hasGeneSelection) {
    return (
      <div className="empty-state">
        Run the gene-selection sweep or the soft power test to see results here.
      </div>
    )
  }

  return (
    <>
      <GeneSelectionPanel
        jobId={geneSelectionJobId}
        rows={geneSelection}
        onUseFraction={onUseFraction}
      />

      {hasSoftPower && (
        <div className="result-block">
          <div className="result-block-header">
            <span>Soft Power Diagnostic</span>
            <span className="result-header-actions">
              <button
                className="btn-sm"
                onClick={() => triggerDownload(softPowerPlotUrl(softPowerJobId), 'soft-power-diagnostic.png')}
              >
                Export PNG
              </button>
              {/* The PDF the pipeline wrote with ggsave(), fetched off disk through the
                  file endpoint — no dedicated route needed. */}
              {outDir && (
                <button
                  className="btn-sm"
                  onClick={() => triggerDownload(fileDownloadUrl(`${outDir}/soft_power_plot.pdf`), 'soft-power.pdf')}
                >
                  Export PDF
                </button>
              )}
            </span>
          </div>
          <div className="result-block-body">
            <img
              className="result-img"
              src={softPowerPlotUrl(softPowerJobId)}
              alt="Soft Power Diagnostic"
              loading="lazy"
            />
          </div>
        </div>
      )}

      <SoftPowerTable softPowers={softPowers} />

      {analyzeJobId && (
        <PlotBlock
          title="kME Plot (Module Hub Genes)"
          src={networkPlotUrl(analyzeJobId)}
          filename="kme-plot.png"
        />
      )}

      {analyzeJobId && hasPreservation && (
        <PlotBlock
          title={`Module Preservation — ${refGroup} (reference) vs ${queryGroup} (query)`}
          src={preservationPlotUrl(analyzeJobId)}
          filename="module-preservation.png"
        />
      )}

      {hasPreservation && (
        <div className="result-block">
          <div className="result-block-header">
            <span>Preservation Statistics — {preservation.length} modules</span>
            <span className="result-header-actions">
              <button className="btn-sm" onClick={() => downloadTable(preservation, 'preservation.csv', ',')}>Export CSV</button>
              <button className="btn-sm" onClick={() => downloadTable(preservation, 'preservation.tsv', '\t')}>Export TSV</button>
            </span>
          </div>
          <div className="result-block-body" style={{ padding: 0 }}>
            <PreservationTable rows={preservation} />
          </div>
          <div className="result-block-footer">
            Sorted least-preserved first. Zsummary &lt; 2 = not preserved in {queryGroup};
            2–10 = weak; &gt; 10 = strongly preserved.
          </div>
        </div>
      )}

      {modules && modules.length > 0 && (
        <div className="result-block">
          <div className="result-block-header">
            <span>
              Module Table — {modules.length} genes
              {hasPreservation && ` (${refGroup} network)`}
            </span>
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
