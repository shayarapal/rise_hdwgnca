import { useState, useMemo } from 'react'

export default function ModuleTable({ rows }) {
  const [sort, setSort] = useState({ col: 'gene_name', dir: 'asc' })

  if (!rows || rows.length === 0) return <p className="empty-state">No module data.</p>

  // Detect columns — fixed first two, then kME_ columns alphabetically
  const allCols = Object.keys(rows[0])
  const kmeCol  = allCols.filter((c) => c.startsWith('kME_')).sort()
  const columns = ['gene_name', 'module_color', ...kmeCol]

  const sorted = useMemo(() => {
    const { col, dir } = sort
    return [...rows].sort((a, b) => {
      const av = a[col] ?? ''
      const bv = b[col] ?? ''
      const cmp = typeof av === 'number' ? av - bv : String(av).localeCompare(String(bv))
      return dir === 'asc' ? cmp : -cmp
    })
  }, [rows, sort])

  const toggle = (col) =>
    setSort((s) => ({ col, dir: s.col === col && s.dir === 'asc' ? 'desc' : 'asc' }))

  const fmt = (val) =>
    typeof val === 'number' ? val.toFixed(3) : val ?? '—'

  return (
    <div className="module-table-wrap">
      <table className="module-table">
        <thead>
          <tr>
            {columns.map((col) => (
              <th key={col} onClick={() => toggle(col)}>
                {col === 'gene_name'    ? 'Gene'   :
                 col === 'module_color' ? 'Module' :
                 col.replace('kME_', 'kME ')}
                {sort.col === col && (
                  <span className="sort-indicator">{sort.dir === 'asc' ? '▲' : '▼'}</span>
                )}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {sorted.map((row, i) => (
            <tr key={i}>
              {columns.map((col) => (
                <td key={col}>
                  {col === 'module_color' ? (
                    <span className="module-chip">
                      <span
                        className="module-chip-dot"
                        style={{ background: row.module_color }}
                        title={row.module_color}
                      />
                      {row.module_color}
                    </span>
                  ) : (
                    fmt(row[col])
                  )}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}
