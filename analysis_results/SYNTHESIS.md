# Synthesis: What This Body of Work Shows, and How to Use It

Nine analyses across two species and two disease models, all built around one question: **is
CACNA1D's co-expression network disrupted in Parkinson's disease, in a way that points to
additional drug targets beyond isradipine (which failed STEADY-PD III)?**

This document is the map: what each study actually shows, its strength of evidence, and where it
belongs in a paper. Read the individual `README.md` in each subfolder for full detail — this is
the cross-study synthesis.

## The nine analyses at a glance

| # | Study | Comparison | Evidence type | CACNA1D's module | Statistically significant? |
|---|---|---|---|---|---|
| 1 | Human GSE243639, SNc | PD vs control (2 separate networks) | Qualitative (independently-built networks) | turquoise (control, kME=0.15) → blue (PD, kME=0.54) | Not a formal test — see #2 |
| 2 | Human GSE243639, SNc | PD vs control (1 combined network) | **Formal DME test** | turquoise, kME=0.53 | **No** (p.adj=1.0) — but 2 *other* modules are (brown p=2.5e-22, blue p=5.0e-4) |
| 3 | Mouse GSE233866, healthy | SNc vs VTA (2 separate networks) | Qualitative | turquoise (SNc, kME=0.44) → brown (VTA, kME=0.30) | Not a formal test — see #4 |
| 4 | Mouse GSE233866, healthy | SNc vs VTA (formal preservation test) | **Formal statistical test** | turquoise, **Zsummary=20.4 (strongly preserved)** | Yes — module structure holds across regions |
| 5 | Mouse GSE233866, 6-OHDA | SNc lesioned vs intact (2 separate networks) | Qualitative | blue (intact, kME=0.35) → brown (lesioned, kME=0.40) | Not a formal test — see #6 |
| 6 | Mouse GSE233866, 6-OHDA | SNc lesioned vs intact (1 combined network) | **Formal DME test** | blue, kME=0.40 | **Yes — extremely** (p.adj=6.1×10⁻²⁶⁹) |
| 7 | Mouse GSE233866, 6-OHDA | SNc vs VTA, **within lesioned cells only** (1 combined network) | **Formal DME test** | blue, kME=0.32 | **Yes — extremely, and highly specific** (p.adj=4.5×10⁻¹²³, higher in SNc) |
| 8 | Mouse GSE233866, 6-OHDA | VTA lesioned vs intact (1 combined network, added 2026-08-13) | **Formal DME test** | green, kME=0.27 | **Yes, but modestly** (p.adj=2.2×10⁻⁵ — same direction as #6, ~10²⁶⁴-fold less extreme) |
| 9 | Mouse GSE233866, 6-OHDA | SNc vs VTA, **within intact cells only** (1 combined network, added 2026-08-13) | **Formal DME test** | turquoise, kME=0.33 | **Yes — extremely, essentially identical to #7** (p.adj≈0, +13.21 log2FC vs. #7's +13.34) |

Study #7 answers a different question than #3/#4 or #5/#6: not baseline regional structure, and
not disease-vs-healthy within one region, but **what's more active in SNc than VTA while
lesioned cells are being compared** — originally framed as the signal most directly relevant to
a region-selective treatment target, specific to neurodegeneration itself. 11 of 12 modules
differ significantly here (expected — SNc and VTA are fundamentally different populations), and
CACNA1D's own module is among the strongest signals in the whole table.

**Study #9 revises that framing, and this is important: the SNc-vs-VTA gap is not
disease-induced.** #9 runs the identical design on the same six animals' *intact* hemispheres
and finds essentially the same result as #7 (+13.21 log2FC vs. +13.34, both p.adj effectively
0) — meaning SNc and VTA are already this molecularly distinct for CACNA1D's program **before**
any lesioning, and lesioning does not detectably widen the gap. Report #7/#9 together as "SNc
and VTA are constitutively distinct for this program, and that distinction persists through
disease" — not "this activates during neurodegeneration." This is arguably a *stronger* basis
for a region-selective target (an anatomical baseline difference, not something contingent on
catching an active disease window), but it is a different claim than what #7 alone would
suggest, and should not be presented as if it were.

**Study #8 is the necessary check on #6's framing, and the honest result changes that framing.**
#6 alone could be (mis)read as "CACNA1D's network response to lesioning is an SNc-specific
phenomenon." #8 tests that directly by running the identical lesioned-vs-intact DME design within
VTA instead, and finds the **same direction** of effect (CACNA1D's module downregulated after
lesioning) at **far lower magnitude and significance** (p.adj=2.2×10⁻⁵ vs. 6.1×10⁻²⁶⁹; ranked 10th
of 12 modules in VTA's table vs. 2nd of 7 in SNc's). The standalone-network kME trajectory shows
the same pattern: VTA's kME also rises with lesioning (0.230→0.283, +23% relative — actually a
slightly larger relative increase than SNc's own +15%). **Do not claim VTA is unaffected by
lesioning — it isn't.** The defensible claim is narrower and more precise: the response exists in
both regions, but is dramatically larger in SNc, and #7 (SNc-vs-VTA compared directly within the
lesioned state) remains the strongest, unaffected evidence that SNc's version of this program is
the one that matters most during active neurodegeneration.

## What each pair (qualitative + formal) actually means

Every comparison was run two ways, and the difference matters:

- **Two separately-built networks** (the "qualitative" rows) let each condition define its own
  modules independently, then compare where CACNA1D lands and how strong its connectivity (kME)
  is in each. This is real, useful evidence that CACNA1D's co-expression *neighborhood*
  reorganizes — but module color/identity is arbitrary per run, so it can't tell you whether a
  specific module's *activity level* differs by condition in a statistically rigorous sense.
- **One combined network + FindDMEs** (the "formal" rows) fixes a single module definition across
  both conditions, then tests whether each module's mean activity (hME) differs between them.
  This is the statistically defensible claim — "module X's activity changes with condition,
  p=..." — and is what a reviewer will actually want to see for a differential-activity claim.

Both are legitimate; they're not redundant, they answer different questions ("does the
neighborhood reorganize" vs. "does average activity change"). Use both, but be precise about
which claim each one supports.

## The central, honest finding

**The human PD-vs-control formal DME test does NOT show CACNA1D's own module as significantly
differentially active.** Two *other* human modules do (brown, down in PD, containing SACS and
AKAP9; blue, up in PD, containing UCHL1 — an established PD gene — plus HSP90AA1/AB1 and
NEFL/NEFM). Meanwhile, **the mouse 6-OHDA lesion model shows CACNA1D's own module as extremely,
overwhelmingly differentially active** (p.adj=6.1×10⁻²⁶⁹).

This is not a contradiction — it's a difference in what an acute chemical lesion vs. chronic
naturally-occurring human disease produces. 6-OHDA is a synchronized, severe experimental insult;
it should produce a much larger, easier-to-detect module-level shift than slowly progressive,
heterogeneous, postmortem human disease. Don't frame this as "the mouse result is better" or "the
human result failed" — frame it as: the mouse model gives a strong, clean signal for CACNA1D's
network responding to acute dopaminergic injury; the human data gives a subtler, more clinically
realistic picture where CACNA1D sits in a network alongside — but not identical to — the modules
that most clearly shift with disease.

## Where each result belongs in your paper

**Motivation / Introduction** — not something this pipeline produced, but the framing: CACNA1D
(Cav1.3) L-type calcium channel pacemaking, SNc-selective vulnerability, isradipine's clinical
failure (STEADY-PD III). This work is testing whether a network-based, multi-gene approach can
succeed where single-target Cav1.3 blockade did not.

**Results — regional vulnerability context** — studies #3/#4. Use these to establish that SNc and
VTA are structurally distinct at the network level even without disease (module identity shift +
the strongly-preserved-neighborhood finding), motivating *why* SNc is the focus. This is
background/mechanistic, not a disease finding — the data has no lesion or disease in it.

**Results — human disease signal (primary translational claim)** — studies #1/#2. Lead with the
formal DME result (#2) for rigor: report which modules ARE significantly differentially active
(brown, blue) and their hub genes (SACS, AKAP9, UCHL1, HSP90AA1/AB1, NEFL/NEFM) as the
statistically supported candidate targets. Use #1 (CACNA1D's own kME shift, 0.15→0.54) as
secondary, clearly-labeled-as-exploratory evidence that CACNA1D's specific neighborhood
reorganizes, even though that module's average activity didn't reach significance.

**Results — mouse disease-model validation** — studies #5/#6. Use #6 (the extremely significant
mouse DME) as same-species, well-powered statistical confirmation that CACNA1D's network
responds to dopaminergic injury — a positive control, in effect, showing this kind of network
perturbation is real and detectable when the insult is large enough. Report the effect-size
difference vs. human openly (acute vs. chronic) rather than letting a reader assume they should
match. Pair this with #8 (the VTA-only version of the same test) rather than presenting #6 in
isolation — report plainly that VTA shows the same direction of response, at markedly lower
magnitude, so a reviewer sees you tested the specificity claim yourself rather than assuming it.

**Results — constitutive region-selective signal (the sharpest single finding)** — studies
#7/#9 together, not #7 alone. CACNA1D's module is one of the most significant differences in
the entire SNc-vs-VTA table both during lesioning (#7, p.adj=4.5×10⁻¹²³) and in intact tissue
(#9, p.adj≈0, essentially the same magnitude) — specifically *more active in SNc*, the region
that's selectively vulnerable, regardless of disease state. Lead with this pair alongside #6 as
the strongest quantitative results in the whole project — but (a) note the caveat that 11 of 12
modules differ here, since SNc and VTA are different populations to begin with, so frame
CACNA1D's result as "notably strong even against that broad background," not as a uniquely
surprising isolated hit, and (b) state plainly that #9 shows this is a baseline/constitutive
difference, not something lesioning switches on — that's a more precise, arguably stronger,
claim than "active during neurodegeneration," but it is a different claim and should be worded
that way. This is also the result that #8 leaves untouched and that should carry the
region-selectivity claim — not #6 read on its own.

**Discussion — candidate targets beyond isradipine** — synthesize across all of it: propose
UCHL1, SACS, AKAP9, HSP90AA1/AB1, PDE4D, PRKG1 (this last pair from the earlier informal
PD-network hub list) as a candidate multi-gene network for follow-up, explicitly noting PDE4D
already has approved inhibitor drugs (roflumilast, apremilast) as a repurposing angle. State the
cross-species evidence as *convergent but not equivalent* — direction of effect (network
perturbation under disease) replicates, specific genes and magnitude do not need to match for the
finding to be meaningful.

**Limitations** — state plainly: (1) CACNA1D's own module did not reach significance in the
rigorous human test; (2) module eigengene log2FC values can be numerically unstable near zero —
several rows in the DME tables have implausibly large log2FC that should not be quoted as real
fold-changes, only the p-values are trustworthy there; (3) sample sizes are modest throughout
(human: 15 PD/14 control donors; mouse: 6 animals); (4) SNc/VTA labels in the mouse data are
marker-derived (Sox6/Calb1), not provided by the original study, and ~30% of DA neurons in the
untreated arm were left unlabeled as ambiguous; (5) cross-species findings are convergent, not a
substitute for direct human validation of any proposed target; (6) VTA's response to lesioning
(#8) is directionally the same as SNc's, not absent — the region-selectivity claim rests on
magnitude and on the direct SNc-vs-VTA comparison (#7), not on VTA showing no response at all,
and should be worded that precisely in any writeup.

## File map

- `../archive/human_GSE243639/analysis_results/pd_vs_control_separate_networks/{control_network,pd_network}/` — study #1 (archived 2026-08-13 pending a rebuild from raw data; see `../archive/human_GSE243639/README.md`)
- `../archive/human_GSE243639/analysis_results/pd_vs_control_combined_dme/` — study #2 (archived, same note as above)
- `mouse_GSE233866/healthy_baseline/{snc_network,vta_network}/` — study #3
- `mouse_GSE233866/healthy_baseline/snc_vs_vta_preservation/` — study #4
- `mouse_GSE233866/lesion_model/{snc_lesioned_network,snc_intact_network}/` — study #5
- `mouse_GSE233866/lesion_model/lesioned_vs_intact_combined_dme_snc/` — study #6
- `mouse_GSE233866/lesion_model/vta_lesioned_network/` + `snc_vs_vta_combined_dme_lesioned/` — study #7
- `mouse_GSE233866/lesion_model/vta_intact_network/` + `lesioned_vs_intact_combined_dme_vta/` — study #8 (added 2026-08-13)
- `mouse_GSE233866/lesion_model/snc_vs_vta_combined_dme_intact/` — study #9 (added 2026-08-13)
- `reminders/reminder.md` — running project log of decisions and pivots made along the way
