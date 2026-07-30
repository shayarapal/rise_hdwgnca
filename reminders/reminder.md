# Reminder — pivot to SNc vs. VTA (2026-07-29)

## Status: PD-vs-control (SNc only) — PAUSED, not abandoned

Everything already built for GSE243639 (SNc dopaminergic neurons, 15 PD / 14 control donors) is
complete and lives in `analysis_results/{control_network,pd_network}/`: dendrogram, kME plot,
hME UMAP, module correlogram, UCell module-score UMAP, Enrichr enrichment. Reproducible via
`DATA_UNZIPPED/build_control_network.R` and `build_pd_network.R`.

Headline finding so far: CACNA1D sits in a weak, non-specifically-enriched module in control
(turquoise, kME=0.15) vs. a strong module enriched for synaptic/calcium-channel signaling in PD
(blue, kME=0.54). The formal statistical preservation test (`ModulePreservation`, PD vs control)
never completed — it kept exceeding available RAM even after fixing two real memory bugs and
raising the ceiling to 44GB. See main `analysis_results/README.md` for full detail.

**Do not delete or redo this work.** Control donors can be reintegrated later (e.g. a four-way
SNc-PD / SNc-control / VTA-PD / VTA-control comparison) if it turns out to matter once the VTA
side exists.

## Pivot: the real comparison is SNc vs. VTA, not PD vs. control

Rationale (confirmed sensible): SNc and VTA are both midbrain dopaminergic populations, but SNc
neurons are selectively vulnerable in Parkinson's while VTA neurons are comparatively spared —
one of the central open questions in PD neuroscience. CACNA1D (Cav1.3, L-type calcium channel)
is specifically implicated in that selective vulnerability via calcium-driven pacemaking stress.
Comparing CACNA1D's co-expression network between SNc and VTA addresses the vulnerability
question directly, more so than a PD-vs-control comparison within SNc alone.

**Blocker:** GSE243639 is SNc-only (stated in its own methods) — it cannot support this
comparison. A new dataset containing VTA (ideally VTA + SNc from the same study, or at least a
comparable VTA dataset with the same cell-type-defining markers) is required.

## Species decision: no human VTA scRNA-seq dataset found

User couldn't find a human single-cell VTA dataset. Decision: use a rodent (mouse) dataset
instead, but ONLY one that profiles **both SNc and VTA in the same study** — directly
cross-comparing rat/mouse VTA against the existing human SNc data (GSE243639) would confound
species differences with region differences and was rejected as not valid for a primary claim.

## Candidate dataset found and verified: GSE233866

"Transcriptomic atlas of midbrain dopamine neurons uncovers differential vulnerability in a
Parkinsonism lesion model" — Yaghmaeian Salmani et al., Karolinska Institute
(https://elifesciences.org/articles/89482).

- Mouse, ~70,000 midbrain nuclei, snRNA-seq
- Explicitly distinguishes SNc vs VTA mDA neurons via Sox6 (SNc) / Calb1 (VTA) — near-mutually-
  exclusive expression defines the split in this dataset
- Also includes a 6-OHDA lesion model (untreated vs. lesioned), so this single dataset could
  later support a disease-relevant comparison too, without needing a third dataset
- Supplementary files on GEO: `GSE233866_untreated_counts.csv.gz` (31.5MB — the baseline SNc/VTA
  atlas, what to start with) and `GSE233866_lesion_intact_counts.csv.gz` (282.7MB — the 6-OHDA
  experiment, for later)

## Progress: GSE233866 Seurat object built (2026-07-29)

`DATA_GSE233866/build_seurat_GSE233866.R` ran successfully. Note: `data.table::fread` silently
mis-parses this file's header (which has one fewer field than the data rows, no name for the
gene-symbol column) — it discards the real header and uses the first gene's data row as column
names instead, losing a gene and all cell barcodes. Fixed by reading the header line separately
and passing `header=FALSE, skip=1` to fread. If rebuilding, do not revert to plain
`fread(file, header=TRUE)`.

Results: 8,311 → 8,183 cells post-QC. 4,514 Dopaminergic Neurons (Th/Slc6a3/Ddc/Slc18a2) vs 3,669
non-DA. Within DA neurons, Sox6/Calb1 split: **2,226 SNc, 832 VTA**, 1,456 left unlabeled
(ambiguous — both or neither marker). Saved as `seurat_GSE233866_SNc_VTA.rds` (+ intermediate
`seurat_GSE233866_clustered.rds`, pre-DA/region-labeling checkpoint).

Also added to r-service/Dockerfile: `harmony` (animal batch correction) and `R.utils`
(fread .gz support) — both were missing.

## Done: SNc and VTA networks built and compared (2026-07-29)

Both built via `DATA_GSE233866/build_snc_network.R` / `build_vta_network.R`. Full writeup in
`analysis_results/mouse_GSE233866/README.md`. Headline: CACNA1D lands in a different module in
each region (SNc: turquoise, kME=0.048; VTA: brown, kME=0.094) but both connections are weak —
much weaker than the human PD-vs-control result — likely a sample-size effect (2,226 / 832 cells
vs. tens of thousands in the human data), not necessarily a weaker true biological signal.
SNc resolved into 9 modules vs. VTA's 39 (many small) — also plausibly a sample-size artifact,
worth rechecking with a larger VTA dataset if one turns up.

Both networks were very cheap to build (SNc peaked at 3.87GB RAM, VTA at 2.25GB) — nowhere near
the memory ceiling that blocked the human PD-vs-control preservation test. This means, unlike the
human data, a formal `ModulePreservation` (SNc-as-ref vs VTA-as-query) permutation test is
actually feasible here on this hardware if wanted next.

## Done: formal ModulePreservation test, SNc vs VTA (2026-07-29)

`DATA_GSE233866/run_preservation_mouse.R`, 200 permutations, completed in ~12 min at 3.88GB peak
(no memory issues at all, unlike the human data). Full writeup in
`analysis_results/mouse_GSE233866/README.md`.

Result: CACNA1D's module (turquoise) has Zsummary=20.4 — strongly preserved (>10 threshold),
ranking 3rd of 8 real modules. Refines the earlier qualitative finding: the gene community
CACNA1D belongs to in SNc clearly holds together in VTA too, even though CACNA1D itself got
assigned a different-colored module (brown) when VTA was analyzed independently — likely a shift
in CACNA1D's own relative connectivity within a stable neighborhood, not a dissolving module.

## Done: formal DME tests, human + mouse lesion model (2026-07-29)

Both use ONE combined network across conditions + `FindDMEs`, not the separately-built networks
(different module definitions, not comparable module-for-module). Human (PD vs control, DA
neurons pre-filtered first — 6.07GB vs the 25-34GB the original per-condition builds needed by
processing all cell types): CACNA1D's own module does NOT reach significance (p.adj=1.0), but two
others do (brown: SACS/AKAP9; blue: UCHL1 — an established PD gene — + HSP90AA1/AB1, NEFL/NEFM).
Mouse (SNc lesioned vs intact): CACNA1D's own module IS extremely significant (p.adj=6.1e-269) —
far more profound than human, consistent with an acute lesion producing a bigger signal than
chronic disease. Full writeup: `analysis_results/SYNTHESIS.md` (new top-level cross-study doc).

Also needed: `ggforestplot` (GitHub, NightingaleHealth/ggforestplot) added to Dockerfile —
`PlotDMEsLollipop` needs it and crashes trying to auto-install itself without `devtools`.

## Done: SNc vs VTA, within lesioned state only (2026-07-29)

User's ask: not baseline region difference, not disease-vs-healthy within one region, but what's
active in SNc vs VTA *specifically during* neurodegeneration — the most direct region-selective
treatment signal. Built `vta_lesioned_network` (qualitative) + `snc_vs_vta_combined_dme_lesioned`
(formal DME, SNc-lesioned vs VTA-lesioned pooled, 4,011 cells, 4.12GB peak).

Result: 11 of 12 modules differ significantly (expected — SNc/VTA are different populations to
begin with), but CACNA1D's own module (blue, kME=0.32) is among the strongest signals in the
whole table (p.adj=4.5e-123) and is specifically *more active in SNc* — arguably the single
sharpest result in the whole project for "why SNc specifically."

## Done: reorganized analysis_results for clarity (2026-07-29)

Old flat structure (`control_network/`, `pd_network/`, `human_snc_dme/` at top level;
`mouse_GSE233866/{snc_network,vta_network,snc_dme,...}/` all flat) reorganized into:
`human_GSE243639/{pd_vs_control_separate_networks,pd_vs_control_combined_dme}/` and
`mouse_GSE233866/{healthy_baseline,lesion_model}/...`, with folder names describing what's
inside rather than needing the README to explain it. New top-level `README.md` is a short index;
`SYNTHESIS.md` remains the full cross-study synthesis. See `SYNTHESIS.md`'s file map for the
complete current path list before assuming an old path still exists.

## Next steps (not yet started)

1. Revisit whether the PD-vs-control angle within VTA specifically (using this dataset's
   6-OHDA lesion arm) adds anything beyond what SNc already covers.
2. Consider whether a larger VTA dataset would change the module-count asymmetry noted in
   `mouse_GSE233866/README.md` (VTA fragments into far more, smaller modules than SNc, likely a
   sample-size artifact rather than a real difference in network organization).
