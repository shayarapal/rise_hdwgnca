import { downloadTable, geneSelectionPlotUrl, triggerDownload } from '../api.js'

const fmt = (v) => (v == null ? '—' : typeof v === 'number' ? v.toLocaleString() : v)

// hdWGCNA has no tuner for `fraction`, so there is no "best" to compute — the sweep just shows
// how many genes each fraction keeps. Picking one is a judgement call: too high and real
// co-expression partners are dropped, too low and the network fills with noise.
export default function GeneSelectionPanel({ jobId, rows, onUseFraction }) {
  if (!jobId || !rows?.length) return null

  const plotSrc = geneSelectionPlotUrl(jobId)

  return (
    <div className="result-block">
      <div className="result-block-header">
        <span>Gene Selection Sweep — genes kept per fraction</span>
        <span className="result-header-actions">
          <button className="btn-sm" onClick={() => triggerDownload(plotSrc, 'gene-selection.png')}>Export PNG</button>
          <button className="btn-sm" onClick={() => downloadTable(rows, 'gene_selection.csv', ',')}>Export CSV</button>
          <button className="btn-sm" onClick={() => downloadTable(rows, 'gene_selection.tsv', '\t')}>Export TSV</button>
        </span>
      </div>

      <div className="result-block-body">
        <img className="result-img" src={plotSrc} alt="Genes selected per fraction" loading="lazy" />
      </div>

      <div className="module-table-wrap">
        <table className="module-table">
          <thead>
            <tr>
              <th>Fraction</th>
              <th>Genes selected</th>
              <th>% of all genes</th>
              <th>Notes</th>
              <th />
            </tr>
          </thead>
          <tbody>
            {rows.map((row) => (
              <tr key={row.fraction} className={row.error ? 'sft-fail' : ''}>
                <td>{row.fraction}</td>
                <td>{fmt(row.n_genes)}</td>
                <td>{row.pct_of_all_genes == null ? '—' : `${row.pct_of_all_genes}%`}</td>
                <td className="gs-note">{row.error || row.warning || ''}</td>
                <td>
                  {row.n_genes > 0 && (
                    <button className="btn-sm" onClick={() => onUseFraction(row.fraction)}>
                      Use this
                    </button>
                  )}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="result-block-footer">
        A gene is kept when it is detected in at least <em>fraction</em> of cells — this is not
        the top X% of genes. hdWGCNA warns below 100 genes and errors at 0; the tutorial default
        is 0.05.
      </div>
    </div>
  )
}
