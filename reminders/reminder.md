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

## Next steps (not yet started)

1. Identify a public dataset with VTA dopaminergic neuron snRNA-seq/scRNA-seq (GEO or similar).
   Likely markers for distinguishing VTA vs SNc DA neurons: SOX6 / ALDH1A1 (SNc-enriched) vs.
   OTX2 / CALB1 / VIP (VTA-enriched, per literature — verify against the chosen dataset's own
   marker annotations rather than assuming).
2. Run the same pipeline (build network → dendrogram → kME → correlogram → hME/UCell UMAP →
   Enrichr) exclusively on VTA neurons from that dataset.
3. Compare CACNA1D's module membership, kME, and enrichment between the VTA network and the
   existing SNc network(s).
4. Only after that: revisit whether the PD-vs-control angle (within SNc, and eventually within
   VTA) is still worth pursuing, and whether DME (differential module eigengene) analysis is the
   right tool for it.
