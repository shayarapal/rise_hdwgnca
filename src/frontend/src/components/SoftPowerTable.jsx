import { downloadTable } from '../api.js'

// GetPowerTable()'s column names contain dots (SFT.R.sq, mean.k.), so they must be indexed
// as strings — row['SFT.R.sq'], never row.SFT.R.sq.
const COLS = [
  { key: 'Power',    label: 'Power' },
  { key: 'SFT.R.sq', label: 'Scale-free R²' },
  { key: 'slope',    label: 'Slope' },
  { key: 'mean.k.',  label: 'Mean k' },
]

const fmt = (v) => (typeof v === 'number' ? v.toFixed(3) : v ?? '—')

export default function SoftPowerTable({ softPowers }) {
  if (!softPowers?.table?.length) return null

  const {
    table, recommended_power, smallest_power, differ,
    sft_threshold, min_power, warning,
  } = softPowers

  return (
    <div className="result-block">
      <div className="result-block-header">
        <span>Soft Power Table — scale-free fit</span>
        <span className="result-header-actions">
          <button className="btn-sm" onClick={() => downloadTable(table, 'soft_power_table.csv', ',')}>Export CSV</button>
          <button className="btn-sm" onClick={() => downloadTable(table, 'soft_power_table.tsv', '\t')}>Export TSV</button>
        </span>
      </div>

      <div className="result-block-body" style={{ padding: 0 }}>
        {warning && <div className="banner banner-warn">{warning}</div>}

        {recommended_power != null && (
          <div className="banner banner-ok">
            Recommended soft power <strong>{recommended_power}</strong> — the lowest power with
            R² ≥ {sft_threshold}{' '}
            {differ ? `above hdWGCNA's min_power of ${min_power}.` : '.'}
          </div>
        )}

        {differ && (
          <div className="banner banner-note">
            The smallest power reaching R² ≥ {sft_threshold} is <strong>{smallest_power}</strong>,
            but ConstructNetwork only ever auto-selects powers above {min_power}, so it would
            build at <strong>{recommended_power}</strong>. Prefilled with {recommended_power};
            override it if you want {smallest_power}.
          </div>
        )}

        <div className="module-table-wrap">
          <table className="module-table">
            <thead>
              <tr>{COLS.map((c) => <th key={c.key}>{c.label}</th>)}</tr>
            </thead>
            <tbody>
              {table.map((row) => {
                const passes = row['SFT.R.sq'] >= sft_threshold
                const isRec  = row.Power === recommended_power
                return (
                  <tr key={row.Power} className={[passes ? 'sft-pass' : '', isRec ? 'sft-rec' : ''].filter(Boolean).join(' ')}>
                    {COLS.map((c) => (
                      <td key={c.key}>
                        {fmt(row[c.key])}
                        {c.key === 'Power' && isRec && <span className="sft-chip">recommended</span>}
                      </td>
                    ))}
                  </tr>
                )
              })}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  )
}
