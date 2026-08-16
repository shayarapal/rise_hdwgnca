# Mouse GSE233866 — CACNA1D Co-expression, SNc vs. VTA and Healthy vs. Lesioned

GSE233866 (Yaghmaeian Salmani et al., Karolinska Institute) profiles mouse midbrain dopaminergic
neurons across two arms: an **untreated** baseline cohort (6 animals, no disease) and a
**6-OHDA lesion** cohort (a different set of 6 animals, each contributing both a lesioned and an
intact/contralateral hemisphere). SNc/VTA labels are not provided by the dataset for either arm —
they're derived here via Sox6 (SNc) / Calb1 (VTA) marker scoring, mirroring the source paper's
own method (see `DATA_GSE233866/build_seurat_GSE233866.R` / `build_seurat_lesion_intact.R`).

Three sub-analyses, each answering a different question. See `analysis_results/SYNTHESIS.md`
for the full cross-study picture. **[`cav13_expression_atlas.html`](cav13_expression_atlas.html)**
(added 2026-08-14) is a self-contained interactive 3D viewer — open it directly in a browser —
plotting all 20,132 SNc/VTA-labeled cells across every network below at once (UMAP1/UMAP2 as
x/y, each cell's own *Cacna1d* expression as height), filterable by condition.

1. **`healthy_baseline/`** — SNc vs. VTA, no disease. *Are the regions structurally different at
   baseline?*
2. **`lesion_model/*_network` + `lesioned_vs_intact_combined_dme_snc/`** — lesioned vs. intact,
   within SNc. *Does disease change this region's network?*
3. **`lesion_model/vta_lesioned_network/` + `snc_vs_vta_combined_dme_lesioned/`** — SNc vs. VTA,
   within the lesioned state only. *What's specifically active in the vulnerable region during
   neurodegeneration itself* — the most direct analog of a region-selective treatment target.
4. **`lesion_model/vta_intact_network/` + `lesioned_vs_intact_combined_dme_vta/`** (added
   2026-08-13) — lesioned vs. intact, within VTA. *Does VTA respond to lesioning the same way
   SNc does, or is SNc's response actually unique to SNc?*
5. **`lesion_model/snc_vs_vta_combined_dme_intact/`** (added 2026-08-13) — SNc vs. VTA, within
   the **intact** state. *Is SNc already this different from VTA before any lesioning, or does
   §3's gap only appear once disease starts?* — see §5 below: it's already there at baseline.
6. **`healthy_baseline/snc_vs_vta_combined_dme_healthy/`** (added 2026-08-16) — SNc vs. VTA,
   formal DME test, **zero disease model involved** (the untreated cohort, not the lesion arm).
   *Is the SNc-vs-VTA activity gap present with no disease model at all, closing the loop
   §5 opened?* — yes: avg_log2FC=+14.55, p.adj=4.17×10⁻²¹², the single strongest result of
   all three SNc-vs-VTA DME tests in this project. See `snc_vs_vta_volcanoes/` for all three
   states plotted together.

---

## 1. `healthy_baseline/` — SNc vs. VTA, no disease

Population sizes: **2,226 SNc cells, 832 VTA cells** (pooled across all 6 untreated animals), out
of 4,514 total dopaminergic neurons identified via Th/Slc6a3/Ddc/Slc18a2 scoring. Soft power = 5
for both (independently determined per network, not assumed from the human dataset).

### Result: CACNA1D's module differs by region; connectivity is moderate in both, stronger in SNc

**Correction (2026-08-13):** this section previously reported kME=0.048 (SNc) and kME=0.094
(VTA), each "very weak." Both were misread from `modules.csv` — pulled from the wrong `kME_*`
column rather than each gene's actual assigned-module column. The correct values, re-verified
directly against `modules.csv` in both folders, are below.

| Region | Module | kME (own module) | Module size (genes) | # modules total |
|---|---|---|---|---|
| SNc | turquoise | 0.442 (moderate) | 1,352 | 9 |
| VTA | brown | 0.296 (weak-to-moderate) | 316 | 39 |

CACNA1D's module membership is a real, non-trivial connection in **both** mouse regions, not
noise — and notably, SNc's connectivity (0.44) is stronger than VTA's (0.30), the opposite
ranking of what this section previously (incorrectly) reported.

This baseline SNc value (0.44) also happens to be numerically higher than every kME in the
lesion-model section below (0.28–0.40) — **do not read that as a "connectivity drops once any
disease modeling happens" trend.** The healthy-baseline cohort (this section) and the lesion-arm
cohort (below) are two entirely separate sets of 6 animals from two different study arms; their
networks were built and clustered independently, from different cell counts, with different
module compositions. They are not a before/after series. The only *paired, controlled*
comparison in this project is intact-vs-lesioned **within** the lesion-arm cohort (0.346→0.398,
same animals, two hemispheres each) — that comparison, and the "connectivity increases after
lesioning" claim built on it, is correct and unaffected by this fix. The human PD-vs-control
comparison (kME 0.15 → 0.54) elsewhere in this project was independently re-checked against its
own `modules.csv` files as part of this fix and is correct as stated — this misread was isolated
to the two mouse healthy-baseline files. Treat the SNc-vs-VTA module identity shift
as informative on its own; don't over-read the exact kME magnitudes given the noise inherent to
per-gene estimates at this cell count. Similarly, the stark module-count difference (SNc: 9 vs.
VTA: 39, many under 100 genes) is likely a sample-size artifact (832 vs. 2,226 cells), not a
confirmed biological difference in network organization.

**Enrichment:** SNc's turquoise module hits "Dopaminergic synapse" (KEGG, p.adj=2.8×10⁻⁴)
directly; VTA's brown module hits "calcium ion transmembrane import into cytosol" (GO,
p.adj=7.8×10⁻³), on-topic for a calcium channel gene but a weaker signal.

### Formal preservation test (`snc_vs_vta_preservation/`, SNc as reference, 200 permutations)

Ran to completion easily (peaked at 3.88GB, ~12 minutes) — the memory ceiling that blocked this
same test on the human data was never a factor here. By WGCNA convention (Langfelder et al.
2011): Zsummary > 10 = strongly preserved, 2–10 = weak/moderate, < 2 = not preserved.

| Module | Zsummary (preservation) | Rank (of 8 real modules) |
|---|---|---|
| brown | 53.1 | 1 |
| green | 47.2 | 2 |
| **turquoise (CACNA1D)** | **20.4** | **3** |
| yellow | 14.9 | 4 |
| blue | 14.6 | 5 |
| red | 10.9 | 6 |
| black | 7.9 | 7 |
| pink | 7.6 | 8 |

(`gold`/`grey` are WGCNA's internal random-control/unassigned-gene bins, excluded from ranking.
`moduleSize.obs` in `preservation.csv` is capped at 1000 genes — a standard WGCNA
computational-tractability limit, not each module's true size.)

**CACNA1D's module (turquoise) is strongly preserved (Zsummary=20.4), ranking 3rd of 8.** This
refines the qualitative finding above: even though CACNA1D gets a *different-colored* label when
VTA is analyzed independently (brown), the broader gene community it belongs to in SNc clearly
still holds together within VTA data — the module-identity shift likely reflects CACNA1D's own
relative connectivity shifting, not that neighborhood dissolving between regions.

All 6 animals cleared the `min_cells=100` metacell threshold, though VTA's smallest sample (s450,
102 cells) was close to that floor.

**Files** — `healthy_baseline/snc_network/`, `vta_network/`: `dendrogram.png`, `network_plot.png`
(kME plot), `module_eigengene_umap.png`, `module_scores_umap.png` (UCell),
`module_correlogram.png`/`.pdf`, `module_by_animal_heatmap.png` (mean hME per module per animal,
z-scored — spots animal outliers), `modules.csv`, `enrichr_table.csv` + `enrichr_plots/*.pdf`,
soft-power diagnostics. `healthy_baseline/snc_vs_vta_preservation/`: `preservation.csv`,
`preservation_plot.png`, `donor_counts.csv`, plus the reference (SNc) network's own outputs.

---

## 2. Lesioned vs. intact within SNc — the mouse PD-model comparison

Uses the study's disease arm: 6-OHDA-lesioned vs. contralateral intact hemisphere, same 6
animals (paired design) — the mouse analog of the human PD-vs-control axis.

Sanity check before any network analysis: 6-OHDA selectively kills dopaminergic neurons, and the
data confirms it — DA neuron counts dropped from 19,662 (intact) to 6,243 (lesioned) overall, and
**SNc lost proportionally more than VTA** (SNc: 8,637→2,291, 3.8x; VTA: 4,426→1,720, 2.6x),
matching this paper's whole premise of differential regional vulnerability.

### Two independently-built networks (`snc_intact_network/`, `snc_lesioned_network/`)

| Condition | n cells | Module | kME (own module) | Peak RAM |
|---|---|---|---|---|
| Intact | 8,637 | blue | 0.346 | 9.98GB |
| Lesioned | 2,291 | brown | 0.398 | 3.82GB |

CACNA1D's connectivity increases modestly under lesioning (0.346→0.398) — same *direction* as
the human PD result, though smaller. Enrichment is similar between the two (axon guidance,
synaptic transmission, cell adhesion) — less differentiated than the human data.

### Formal DME test (`lesioned_vs_intact_combined_dme_snc/`)

One combined network (both conditions pooled, 10,928 cells, 7.67GB peak), then `FindDMEs`.
CACNA1D's module here: **blue, kME=0.40**.

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| yellow | +6.93 | ~0 | ↑↑ lesioned |
| **blue (CACNA1D)** | **−1.66** | **6.1×10⁻²⁶⁹** | **↓ lesioned, extremely significant** |
| red | −1.85 | 1.1×10⁻¹⁴¹ | ↓ lesioned |
| brown | (unstable, see caveat) | 7.1×10⁻¹⁷ | ↓ lesioned |
| turquoise | +1.82 | 2.1×10⁻¹⁶ | ↑ lesioned |
| black | −2.79 | 2.2×10⁻⁴ | ↓ lesioned |
| green | +0.75 | 1.0 (n.s.) | no change |

**CACNA1D's module is massively downregulated after lesioning** (p.adj ≈ 6×10⁻²⁶⁹) — far more
statistically extreme than the equivalent human DME result
(`archive/human_GSE243639/analysis_results/pd_vs_control_combined_dme/` — archived
2026-08-13 pending a rebuild, see `archive/human_GSE243639/README.md`; CACNA1D's module does
NOT reach significance there). Combined with
the kME increase above: CACNA1D's broader co-expression program is suppressed in surviving
lesioned neurons, but CACNA1D itself becomes *relatively* more central to whatever remains of
that shrinking program — a hypothesis worth stating, not an established mechanism.

The brown module's log2FC is numerically unstable (module eigengenes can be negative/near-zero) —
trust the p-value, not that magnitude. Same caveat applies throughout every DME table in this
project.

---

## 3. SNc vs. VTA, within the lesioned state — what's active in the vulnerable region *during* disease

**§5 revises this section's framing (added 2026-08-13):** the SNc-vs-VTA gap below turns out to
already be present in the *intact* hemisphere at essentially the same magnitude — it is a
standing baseline difference, not something that emerges during neurodegeneration. Read
"during disease" below as "present at baseline and preserved through disease," not
"disease-induced." The numbers themselves are unaffected.

The most disease-and-region-specific comparison: not baseline regional structure (§1), not
disease-vs-healthy within one region (§2), but which genes/modules are more active in SNc than
VTA specifically while neurodegeneration is happening — the signal most directly relevant to a
region-selective treatment target.

**`lesion_model/vta_lesioned_network/`** (qualitative, mirrors §1/§2's separately-built
networks): CACNA1D lands in module blue, kME=0.283, 1,720 cells, 4.14GB peak.

**`lesion_model/snc_vs_vta_combined_dme_lesioned/`** (formal DME): one combined network pooling
SNc-lesioned (2,291 cells) + VTA-lesioned (1,720 cells) = 4,011 cells, 4.12GB peak. CACNA1D's
module here: **blue, kME=0.32**.

`FindDMEs(barcodes1=SNc, barcodes2=VTA)` — positive log2FC means higher in SNc:

| Module | avg_log2FC | p.adj | Direction |
|---|---|---|---|
| greenyellow | −5.34 | ~0 | higher in VTA |
| turquoise | −6.70 | 7.2×10⁻²⁸⁸ | higher in VTA |
| tan | +4.70 | 1.1×10⁻²⁷¹ | higher in SNc |
| magenta | +6.83 | 3.4×10⁻¹⁴² | higher in SNc |
| **blue (CACNA1D)** | **+13.34** | **4.5×10⁻¹²³** | **higher in SNc, extremely significant** |
| pink | +5.41 | 2.0×10⁻¹¹⁹ | higher in SNc |
| green | +4.45 | 1.0×10⁻⁵¹ | higher in SNc |
| purple | −0.90 | 4.0×10⁻³² | higher in VTA |
| black | −1.36 | 3.4×10⁻³⁰ | higher in VTA |
| yellow | +12.78 | 8.8×10⁻²⁵ | higher in SNc |
| red | +0.28 | 1.5×10⁻¹⁷ | higher in VTA |
| brown | −9.43 | 0.068 (n.s.) | — |

**11 of 12 modules differ significantly between SNc and VTA during lesioning** — expected, since
these are fundamentally different neuron populations, not just disease-state variants of the
same one. What matters is that **CACNA1D's own module is among the strongest signals in the
entire table** (p.adj=4.5×10⁻¹²³) and is specifically *more active in SNc than VTA* during the
neurodegenerative process — the most direct answer this project has produced to "what's driving
SNc's selective vulnerability while it's actually happening." As with every DME table here,
treat the log2FC magnitudes as directional evidence, not literal fold-changes (module eigengenes
can be negative/near-zero, making that math numerically unstable) — the p-values are what to
trust.

---

## 4. VTA's own lesioned vs. intact response — is SNc's story unique to SNc? (added 2026-08-13)

Section 2 shows CACNA1D's module becomes more central (kME 0.346→0.398) and significantly
perturbed at the eigengene level (p.adj=6.1×10⁻²⁶⁹) after lesioning, **within SNc**. Until
this section was added, VTA had no equivalent test — only a standalone lesioned network
(§3), never paired with a VTA-only intact network. That gap matters: without it, "CACNA1D
responds to lesioning" could not be distinguished from "CACNA1D responds to lesioning
*specifically in SNc*."

**`lesion_model/vta_intact_network/`**: 4,426 cells, soft power 4. CACNA1D lands in module
pink, kME=0.230, 132 genes — the smallest, lowest-kME CACNA1D module in this project, and
the only one with zero significant enrichment hits (best p.adj=0.120). Gene overlap with
`vta_lesioned_network`'s module is 14% (19/132 genes), well below SNc's equivalent 48%.

**`lesion_model/lesioned_vs_intact_combined_dme_vta/`** (formal DME, 6,146 pooled cells):
CACNA1D's module (green, 191 genes) is **significantly downregulated after lesioning**
(avg_log2FC=−1.23, p.adj=2.2×10⁻⁵) — the same direction as SNc. The standalone kME also
rises the same direction as SNc: 0.230→0.283 (+23% relative, vs. SNc's +15%).

**Honest conclusion: VTA shows the same qualitative pattern as SNc, not a null result.**
"Only SNc responds to lesioning" is not supported and should not be claimed. What *is*
still true: SNc's effect is far larger — ~10²⁶⁴-fold more extreme by p.adj, numerically
bigger by log2FC, and ranked 2nd-most-extreme of 7 modules in SNc's table vs. 10th of 12 in
VTA's — despite SNc's pooled network having under 2× VTA's cell count, so this gap isn't
just statistical power. And critically, §3's direct SNc-vs-VTA-during-lesioning test
(p.adj=4.5×10⁻¹²³, CACNA1D's module far more active in SNc) is untouched by this and remains
the strongest evidence for region-selectivity. The defensible framing going forward: this
mechanism isn't unique to SNc, but its consequences are — SNc's version of the program is
dramatically larger and more active throughout the lesioned state, which is still a
meaningful basis for a region-selective treatment argument, just a more precise one than
before this section existed.

---

## 5. SNc vs. VTA, within the intact state — is §3's gap disease-induced or already there? (added 2026-08-13)

§3 found CACNA1D's module dramatically more active in SNc than VTA **during lesioning**
(p.adj=4.5×10⁻¹²³), framed as evidence of something specific to the vulnerable region during
neurodegeneration. This section tests the missing counterfactual on the same 6 animals'
**intact** hemispheres: is SNc already this different from VTA before any lesion, or does the
gap only appear once disease starts?

**`lesion_model/snc_vs_vta_combined_dme_intact/`** (formal DME, 13,063 pooled cells: 8,637 SNc
+ 4,426 VTA): CACNA1D's module (turquoise, 1,629 genes, kME=0.327) is higher in SNc with
avg_log2FC=+13.21, p.adj≈0 (underflows double precision) — essentially identical in direction
and magnitude to §3's lesioned result (+13.34, p.adj=4.5×10⁻¹²³).

**Revised conclusion: the SNc-vs-VTA gap is not disease-induced — it's a standing baseline
feature of this cohort that persists through lesioning, not something that emerges during
neurodegeneration.** This changes §3's framing: describe it as "SNc and VTA are already
molecularly distinct for this program before disease, and that distinction remains through
lesioning," not "this activates during neurodegeneration." This doesn't weaken the
region-selective treatment argument — a constitutive anatomical distinction is arguably a more
reliable basis for selectivity than a disease-triggered one — but it changes what the evidence
supports. The genuinely lesioning-*specific* findings remain §2/§4 (does each region's own
network change with lesioning — yes for SNc strongly, yes for VTA weakly), not §3/§5 (which
show a real but constitutive regional difference, present with or without lesioning).

Turquoise here shares 74% of its genes with `snc_intact_network` and 52% with the
healthy-baseline SNc network — the fifth independently-built network in this project to
recover substantially the same SNc gene neighborhood.

---

## 6. SNc vs. VTA, within healthy baseline — closing the loop with zero disease model (added 2026-08-16)

§5 answered "is the gap already there in *intact* lesion-arm tissue" — yes. This section asks
the sharper version: is it there with **no lesion-arm cohort involved at all**, using the
completely separate untreated cohort instead (`healthy_baseline/`, 6 different animals, no
6-OHDA anywhere in the pipeline)?

**`healthy_baseline/snc_vs_vta_combined_dme_healthy/`** (formal DME, 3,058 pooled cells: 2,226
SNc + 832 VTA): CACNA1D's module (yellow) is higher in SNc with avg_log2FC=+14.55,
p.adj=4.17×10⁻²¹² — the single most significant module in an 11-module table, and the
largest-magnitude of all three SNc-vs-VTA DME results in this project (healthy: +14.55,
intact: +13.21, lesioned: +13.34).

**This closes the loop §5 opened.** Three independently-built, formally-tested networks — one
with no disease model whatsoever — all agree: CACNA1D's module is significantly more active in
SNc than VTA, at comparable magnitude, regardless of disease state. The region-selectivity
claim for this project should be stated as a constitutive property of these two regions, full
stop — not qualified by "during lesioning" or "in the lesion-arm cohort." See
`snc_vs_vta_volcanoes/` for all three results plotted together with matched axes.

---

## Caveats (apply throughout this folder)

- Sox6/Calb1-only cells were left unlabeled rather than forced into a region (1,456 of 4,514 DA
  neurons in the untreated arm) — a stricter or looser labeling rule would shift population sizes.
- Enrichr's `GO_Biological_Process_2023` is queried directly with mouse gene symbols (GO isn't
  species-forked in Enrichr); `KEGG_2019_Mouse` (fixed 2026-08-13, was previously the
  wrong-species `KEGG_2021_Human`) is the correct mouse-specific library.
- These are mouse findings. Translating them back to the human results in
  `archive/human_GSE243639/analysis_results/` (archived 2026-08-13 pending a rebuild) is a
  separate, later inference — not something to claim directly from any comparison in this folder.
- Every kME value in this folder was independently re-verified against `modules.csv` on
  2026-08-13 after two healthy-baseline figures were found misread (see §1's correction note).
  If you're citing a kME from an older copy of this README, re-check it against the source file.
