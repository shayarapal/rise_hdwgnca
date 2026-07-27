// Zsummary convention (Langfelder et al. 2011): > 10 strong evidence of preservation,
// 2–10 weak/moderate, < 2 no evidence. A LOW Zsummary is the interesting result — it
// means the module's co-expression structure does not hold in the query condition.
function zClass(z) {
  if (z == null || isNaN(Number(z))) return ''
  const n = Number(z)
  if (n < 2)  return 'z-none'
  if (n < 10) return 'z-weak'
  return 'z-strong'
}

const Z_COL = 'Zsummary.pres'

export default function PreservationTable({ rows }) {
  if (!rows || rows.length === 0) return null

  // Column names come from WGCNA::modulePreservation via GetModulePreservation(); render
  // whatever it returned rather than hardcoding a schema, and only special-case Zsummary.
  const cols = Object.keys(rows[0])

  return (
    <table className="module-table">
      <thead>
        <tr>{cols.map((c) => <th key={c}>{c}</th>)}</tr>
      </thead>
      <tbody>
        {rows.map((row, i) => (
          <tr key={row.module ?? i}>
            {cols.map((c) => (
              <td key={c} className={c === Z_COL ? zClass(row[c]) : ''}>
                {typeof row[c] === 'number' ? row[c].toFixed(3) : String(row[c])}
              </td>
            ))}
          </tr>
        ))}
      </tbody>
    </table>
  )
}
