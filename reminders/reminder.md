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

## Next steps (not yet started)

1. Download and stage GSE233866 (`untreated_counts.csv.gz` first).
2. Build Seurat object, using this dataset's own Sox6/Calb1-based (or however it's annotated)
   SNc/VTA split rather than assuming a marker set from a different study.
3. Run the same pipeline (build network → dendrogram → kME → correlogram → hME/UCell UMAP →
   Enrichr) separately on the SNc and VTA populations within this mouse dataset.
4. Compare CACNA1D's module membership, kME, and enrichment between mouse SNc and mouse VTA.
   Note as a caveat throughout: this is a same-species (mouse) finding — translating it back to
   the human SNc results already in `analysis_results/` is a separate, later inference, not
   something to claim directly.
5. Only after that: revisit whether the PD-vs-control angle (within SNc, and eventually within
   VTA, now possible via this dataset's own 6-OHDA lesion arm) is still worth pursuing, and
   whether DME (differential module eigengene) analysis is the right tool for it.
