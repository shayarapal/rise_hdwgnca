# UI Design Notes — hdWGCNA Bridge Frontend

The React UI should mirror the feel of the hdWGCNA tutorial workflow: linear, stage-by-stage,
with a single break point where the user must make a biological decision before proceeding.

---

## Overall layout

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  hdWGCNA Bridge                                                [job history] │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                               │
│  ── Pipeline stepper (horizontal, top of page) ──────────────────────────── │
│  [1 Setup] → [2 Metacells] → [3 SetDatExpr] → [4 Soft Power ●] →           │
│                                                  [5 Network] → [6 Modules]  │
│                                                                               │
│  ┌─────────────────────────┐   ┌──────────────────────────────────────────┐ │
│  │  Parameters (left panel)│   │  Results (right panel)                   │ │
│  │                         │   │                                           │ │
│  │  h5seurat_path: [____]  │   │  [Soft power diagnostic PNG]             │ │
│  │  out_dir:       [____]  │   │  [Network plot PNG]                       │ │
│  │  cell_type_col: [____]  │   │  [Module table — sortable grid]          │ │
│  │  group_by:      [____]  │   │                                           │ │
│  │  group_name:    [____]  │   │                                           │ │
│  │  wgcna_name:    [____]  │   │                                           │ │
│  │  network_type:  [signed]│   │                                           │ │
│  │                         │   │                                           │ │
│  │  [Run Soft Power Test]  │   │                                           │ │
│  │                         │   │                                           │ │
│  │  ── Soft power break ── │   │                                           │ │
│  │  (unlocked after plot)  │   │                                           │ │
│  │  soft_power: [__]       │   │                                           │ │
│  │                         │   │                                           │ │
│  │  [Run Full Pipeline]    │   │                                           │ │
│  └─────────────────────────┘   └──────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Pipeline stepper

- 6 steps, horizontal, connected with arrows
- Active step highlighted; completed steps have a checkmark
- Step 4 (Soft Power) has a distinct "break point" indicator (e.g., orange dot)
- Steps 5 and 6 are visually locked (greyed) until `soft_power` is entered

---

## Break point — step 4

After "Run Soft Power Test" is clicked:
1. The soft power diagnostic PNG appears in the results panel
2. A numeric input for `soft_power` appears in the parameters panel below the PNG
3. The stepper advances to step 4 (active)
4. Steps 5–6 remain locked until a value is entered
5. Once the user enters `soft_power` and clicks "Run Full Pipeline", steps 5–6 unlock

This mirrors the actual workflow: the researcher inspects the scale-free topology fit
curve and picks the lowest power at which R² ≥ 0.8 (or similar criterion).

---

## Results panel

| Result | Shown after | Format |
|---|---|---|
| Soft power plot | `/test-soft-powers` job done | `<img>` tag pointing at `/results/{job_id}/soft-power-plot` |
| Network plot | `/analyze` job done | `<img>` tag pointing at `/results/{job_id}/plot` |
| Module table | `/analyze` job done | Sortable HTML table or equivalent; columns: `gene_name`, `module_color`, `module_kME`, `hub_gene_score` |

---

## Job status indicator

While a job is running, show a per-stage progress tracker:

```
SetupForWGCNA       ✓ done
MetacellsByGroups   ⟳ running   ← current stage (spinner)
NormalizeMetacells  ○ pending
SetDatExpr          ○ pending
TestSoftPowers      ○ pending
```

This requires the r-service to emit stage-level progress. Currently `/status/{job_id}`
returns only `running | done | failed`. To support stage-level progress, add a
`stage` field to the job store and update it during pipeline execution.

If stage-level progress is not implemented, a simple "Running… / Done / Failed"
indicator is acceptable for v0.1.

---

## Module table columns

Mirror PyWGCNA output schema so the UI feels familiar to WGCNA Python users:

| Column (CSV) | Display label | Type |
|---|---|---|
| `gene_name` | Gene | string |
| `module_color` | Module | coloured chip (hex background) |
| `kME_{module}` | kME | numeric, 3 decimal places |
| `hub_gene_score` | Hub score | numeric; only present after `GetHubGenes` if added |

---

## Technology constraints

- Vite + React 18, no router (single page), no state management library
- No component library unless explicitly requested — plain HTML + CSS is fine for v0.1
- Vite dev proxy configured in `vite.config.js` to forward API calls to `localhost:8200`
- `fetch` (native) for API calls; no Axios
