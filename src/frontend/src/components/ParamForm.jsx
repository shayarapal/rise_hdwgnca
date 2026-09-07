import { useState } from 'react'
import FileBrowser from './FileBrowser.jsx'

export default function ParamForm({ params, onChange, onRunGeneSelection, geneSelectionRunning }) {
  const set = (field) => (e) => onChange({ ...params, [field]: e.target.value })
  const setNum = (field) => (e) => onChange({ ...params, [field]: Number(e.target.value) })
  const setGroupBy = (idx) => (e) => {
    const updated = [...params.group_by]
    updated[idx] = e.target.value
    onChange({ ...params, group_by: updated })
  }

  // Which field the file browser is currently picking for: null | 'h5seurat_path' | 'out_dir'
  const [browsing, setBrowsing] = useState(null)

  // "custom" is a legal SelectNetworkGenes value but the API rejects it (no gene_list param),
  // so it is not offered here.
  const geneSelect = params.gene_select ?? 'fraction'

  // The sweep loads the object and reads its metadata, so it needs the same fields the
  // pipeline does — minus anything about the network itself.
  const canSweep = Boolean(
    params.h5seurat_path?.trim() &&
    params.out_dir?.trim() &&
    params.cell_type_col?.trim() &&
    !geneSelectionRunning,
  )

  // Condition-comparison mode is derived from the params themselves rather than held
  // as separate state, so clearing it also clears the fields the backend keys off.
  const compare = params.condition_col != null
  const toggleCompare = (e) =>
    onChange(
      e.target.checked
        ? { ...params, condition_col: 'condition', ref_group: '', query_group: '' }
        : { ...params, condition_col: null, ref_group: null, query_group: null },
    )

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

      <p className="form-section-title">Compare conditions</p>

      <div className="form-group">
        <label className="checkbox-label">
          <input
            type="checkbox"
            checked={compare}
            onChange={toggleCompare}
          />
          Compare two conditions (module preservation)
        </label>
        <span className="hint">
          Builds the network in the reference condition, projects it into the query
          condition, and scores whether each module survives
        </span>
      </div>

      {compare && (
        <>
          <div className="form-group">
            <label htmlFor="condition_col">Condition column</label>
            <input
              id="condition_col"
              value={params.condition_col ?? ''}
              onChange={set('condition_col')}
              placeholder="condition"
            />
            <span className="hint">Metadata column holding the two groups</span>
          </div>

          <div className="form-row">
            <div className="form-group">
              <label htmlFor="ref_group">Reference</label>
              <input
                id="ref_group"
                value={params.ref_group ?? ''}
                onChange={set('ref_group')}
                placeholder="control"
              />
            </div>
            <div className="form-group">
              <label htmlFor="query_group">Query</label>
              <input
                id="query_group"
                value={params.query_group ?? ''}
                onChange={set('query_group')}
                placeholder="PD"
              />
            </div>
          </div>
          <span className="hint">
            The network is built in the reference and tested in the query
          </span>
        </>
      )}

      <p className="form-section-title">Gene selection</p>

      <div className="form-group">
        <label htmlFor="gene_select">Gene select</label>
        {/* Read with ?? so an untouched form omits the key entirely and the r-service applies
            its own default — the request body stays identical to what it sends today. */}
        <select id="gene_select" value={geneSelect} onChange={set('gene_select')}>
          <option value="fraction">fraction</option>
          <option value="variable">variable</option>
          <option value="all">all</option>
        </select>
        <span className="hint">
          {geneSelect === 'fraction'
            ? 'Genes detected in at least the given share of cells'
            : geneSelect === 'variable'
              ? 'VariableFeatures() — variance across ALL cell types, not within this one'
              : 'Every gene in the assay'}
        </span>
      </div>

      {geneSelect === 'fraction' && (
        <div className="form-group">
          <label htmlFor="fraction">Fraction</label>
          <input
            id="fraction"
            type="number"
            step={0.01}
            min={0.01}
            max={1}
            value={params.fraction ?? 0.05}
            onChange={setNum('fraction')}
          />
          <span className="hint">
            Keeps genes expressed in at least this share of cells — not the top X% of genes.
            Tutorial default is 0.05.
          </span>
        </div>
      )}

      <button
        type="button"
        className="btn btn-secondary"
        onClick={onRunGeneSelection}
        disabled={!canSweep}
      >
        {geneSelectionRunning ? 'Sweeping…' : 'Run gene-selection sweep'}
      </button>
      <span className="hint">
        Counts how many genes each candidate fraction keeps, so you can pick one from data
        rather than by feel. hdWGCNA ships no tuner for this.
      </span>

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
