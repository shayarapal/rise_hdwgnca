import { useState } from 'react'
import FileBrowser from './FileBrowser.jsx'

export default function ParamForm({ params, onChange }) {
  const set = (field) => (e) => onChange({ ...params, [field]: e.target.value })
  const setNum = (field) => (e) => onChange({ ...params, [field]: Number(e.target.value) })
  const setGroupBy = (idx) => (e) => {
    const updated = [...params.group_by]
    updated[idx] = e.target.value
    onChange({ ...params, group_by: updated })
  }

  // Which field the file browser is currently picking for: null | 'h5seurat_path' | 'out_dir'
  const [browsing, setBrowsing] = useState(null)

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>
      <p className="form-section-title">Data paths</p>

      <div className="form-group">
        <label htmlFor="h5seurat_path">Seurat path (.rds or .h5Seurat)</label>
        <div className="input-with-btn">
          <input
            id="h5seurat_path"
            value={params.h5seurat_path}
            onChange={set('h5seurat_path')}
            placeholder="/shared/data/brain.rds"
          />
          <button type="button" className="btn-browse" onClick={() => setBrowsing('h5seurat_path')}>
            Browse…
          </button>
        </div>
        <span className="hint">Pick a file from your computer, or type a path</span>
      </div>

      <div className="form-group">
        <label htmlFor="out_dir">Output directory</label>
        <div className="input-with-btn">
          <input
            id="out_dir"
            value={params.out_dir}
            onChange={set('out_dir')}
            placeholder="/shared/results/run1"
          />
          <button type="button" className="btn-browse" onClick={() => setBrowsing('out_dir')}>
            Browse…
          </button>
        </div>
        <span className="hint">Created automatically if it doesn't exist</span>
      </div>

      {browsing && (
        <FileBrowser
          mode={browsing === 'out_dir' ? 'dir' : 'file'}
          title={browsing === 'out_dir' ? 'Choose an output folder' : 'Choose an h5Seurat file'}
          onClose={() => setBrowsing(null)}
          onSelect={(path) => { onChange({ ...params, [browsing]: path }); setBrowsing(null) }}
        />
      )}

      <p className="form-section-title">Cell grouping</p>

      <div className="form-group">
        <label htmlFor="cell_type_col">Cell type column</label>
        <input
          id="cell_type_col"
          value={params.cell_type_col}
          onChange={set('cell_type_col')}
          placeholder="cell_type"
        />
        <span className="hint">Metadata column in Seurat object</span>
      </div>

      <div className="form-row">
        <div className="form-group">
          <label htmlFor="group_by_0">Group by [0]</label>
          <input
            id="group_by_0"
            value={params.group_by[0] ?? ''}
            onChange={setGroupBy(0)}
            placeholder="cell_type"
          />
        </div>
        <div className="form-group">
          <label htmlFor="group_by_1">Group by [1]</label>
          <input
            id="group_by_1"
            value={params.group_by[1] ?? ''}
            onChange={setGroupBy(1)}
            placeholder="sample_id"
          />
        </div>
      </div>

      <div className="form-group">
        <label htmlFor="group_name">Group name</label>
        <input
          id="group_name"
          value={params.group_name}
          onChange={set('group_name')}
          placeholder="Excitatory Neurons"
        />
        <span className="hint">Cell type to build the co-expression network for</span>
      </div>

      <p className="form-section-title">Network parameters</p>

      <div className="form-group">
        <label htmlFor="wgcna_name">WGCNA name</label>
        <input
          id="wgcna_name"
          value={params.wgcna_name}
          onChange={set('wgcna_name')}
          placeholder="tutorial"
        />
      </div>

      <div className="form-group">
        <label htmlFor="network_type">Network type</label>
        <select id="network_type" value={params.network_type} onChange={set('network_type')}>
          <option value="signed">signed</option>
          <option value="unsigned">unsigned</option>
          <option value="signed hybrid">signed hybrid</option>
        </select>
      </div>

      <div className="form-row">
        <div className="form-group">
          <label htmlFor="k">k (metacells)</label>
          <input id="k" type="number" min={5} max={100} value={params.k} onChange={setNum('k')} />
        </div>
        <div className="form-group">
          <label htmlFor="max_shared">max shared</label>
          <input id="max_shared" type="number" min={1} max={50} value={params.max_shared} onChange={setNum('max_shared')} />
        </div>
      </div>
    </div>
  )
}
