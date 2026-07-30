# Human GSE243639 — CACNA1D Co-expression, PD vs. Control (SNc Dopaminergic Neurons)

Data: GSE243639 (29 donors: 15 Parkinson's, 14 control; substantia nigra pars compacta snRNA-seq).
Pipeline: `r-service` hdWGCNA bridge (`SetupForWGCNA` → `MetacellsByGroups` → `NormalizeMetacells` →
`SetDatExpr` → `TestSoftPowers` → `ConstructNetwork` → `ModuleEigengenes` → `ModuleConnectivity`),
restricted to the "Dopaminergic Neurons" cell type. Soft power = 5 throughout (network_type=
"signed", gene_select="fraction" 0.05, k=25, max_shared=15).

Two sub-analyses here, answering different questions — see `analysis_results/SYNTHESIS.md` for
how they fit together:

## `pd_vs_control_separate_networks/` — qualitative, independently-built networks

Each condition gets its own network (different module identities/colors, not directly
comparable module-for-module). Shows where CACNA1D lands and its connectivity in each.

| Condition | Module   | kME (own module) | Module size (genes) |
|-----------|----------|-------------------|----------------------|
| Control   | turquoise | 0.147 (weak)      | 3,885                |
| PD        | blue      | 0.540 (strong)    | 1,511                |

CACNA1D is only weakly tied to its assigned module in control tissue, but is a strongly-connected
member of the "blue" module in PD tissue — a qualitative shift in its co-expression neighborhood
between conditions.

**Enrichment reinforces this:** CACNA1D's PD-condition module (blue, kME=0.54) is enriched almost
entirely for **synaptic signaling and calcium channel regulation** — directly on-topic for a
voltage-gated calcium channel gene:

| Term | Adjusted p-value |
|---|---|
| Regulation Of Trans-Synaptic Signaling | 6.1×10⁻⁴ |
| Modulation Of Chemical Synaptic Transmission | 9.0×10⁻⁴ |
| Regulation Of Cation Channel Activity | 1.7×10⁻³ |
| **Regulation Of Voltage-Gated Calcium Channel Activity** | 1.7×10⁻³ |
| Regulation Of Calcium Ion Transmembrane Transporter Activity | 2.0×10⁻³ |
| Regulation Of Synapse Assembly | 9.9×10⁻³ |

By contrast, CACNA1D's control-condition module (turquoise, kME=0.15) shows only weak,
non-specific enrichment (actin filament organization, lipid biosynthesis — all tied at the same
borderline p.adj≈0.028).

Files per condition (`control_network/`, `pd_network/`): `modules.csv`, `dendrogram.png`,
`network_plot_kME.png`, `module_eigengene_umap.png` (hME — use for statistical analyses, it's
continuous and comparable across cells), `module_scores_umap.png` (UCell — use for
visualization, more outlier-robust), `module_correlogram.png`/`.pdf`, `enrichr_table.csv` +
`enrichr_plots/*.pdf`, soft-power diagnostics.

## `pd_vs_control_combined_dme/` — the formal, rigorous test

One combined network (both conditions pooled, DA neurons pre-filtered — 10,388 cells, 6.07GB
peak), then `FindDMEs` tests whether each module's *average activity* differs by condition
within that single fixed module definition — the statistically correct way to test "does this
module's activity change," unlike the qualitative comparison above.

**Result: CACNA1D's own module (turquoise, kME=0.53) does NOT reach significance (p.adj=1.0).**
Two *other* modules do: brown (down in PD, p=2.5×10⁻²², contains SACS/AKAP9) and blue (up in PD,
p=5.0×10⁻⁴, contains **UCHL1** — an established PD gene — plus HSP90AA1/AB1, NEFL/NEFM). See
`analysis_results/SYNTHESIS.md` for the full interpretation and how this reconciles with the
qualitative finding above (short version: they answer different questions — "does the
neighborhood reorganize" vs. "does average activity change" — and aren't in conflict).

Files: `modules.csv`, `dendrogram.png`, `network_plot.png`, `dme_results.csv`, `dme_volcano.png`,
`dme_lollipop.png`, soft-power diagnostics.

## Note on memory

Both analyses run well within resource limits once DA neurons are pre-filtered before network
construction. An earlier attempt at a formal `ModulePreservation()` permutation test (a
different, heavier statistical test than `FindDMEs`) was abandoned after repeatedly exceeding
44GB RAM — two real memory bugs were found and fixed along the way (a hdWGCNA sparse-matrix
inefficiency in `SelectNetworkGenes`, patched in `r-service/seuratdisk_compat.R`; and a
`plumber.R` inefficiency holding both conditions' raw data in memory unnecessarily early), but
what remained was a genuine memory requirement of that specific permutation test on this dataset
size. `FindDMEs` (used above) answers a closely related question without that cost.
