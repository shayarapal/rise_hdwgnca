# Synthesis: What This Body of Work Shows, and How to Use It

Seven analyses across two species and two disease models, all built around one question: **is
CACNA1D's co-expression network disrupted in Parkinson's disease, in a way that points to
additional drug targets beyond isradipine (which failed STEADY-PD III)?**

This document is the map: what each study actually shows, its strength of evidence, and where it
belongs in a paper. Read the individual `README.md` in each subfolder for full detail — this is
the cross-study synthesis.

## The seven analyses at a glance

| # | Study | Comparison | Evidence type | CACNA1D's module | Statistically significant? |
|---|---|---|---|---|---|
| 1 | Human GSE243639, SNc | PD vs control (2 separate networks) | Qualitative (independently-built networks) | turquoise (control, kME=0.15) → blue (PD, kME=0.54) | Not a formal test — see #2 |
| 2 | Human GSE243639, SNc | PD vs control (1 combined network) | **Formal DME test** | turquoise, kME=0.53 | **No** (p.adj=1.0) — but 2 *other* modules are (brown p=2.5e-22, blue p=5.0e-4) |
| 3 | Mouse GSE233866, healthy | SNc vs VTA (2 separate networks) | Qualitative | turquoise (SNc, kME=0.05) → brown (VTA, kME=0.09) | Not a formal test — see #4 |
| 4 | Mouse GSE233866, healthy | SNc vs VTA (formal preservation test) | **Formal statistical test** | turquoise, **Zsummary=20.4 (strongly preserved)** | Yes — module structure holds across regions |
| 5 | Mouse GSE233866, 6-OHDA | SNc lesioned vs intact (2 separate networks) | Qualitative | blue (intact, kME=0.35) → brown (lesioned, kME=0.40) | Not a formal test — see #6 |
| 6 | Mouse GSE233866, 6-OHDA | SNc lesioned vs intact (1 combined network) | **Formal DME test** | blue, kME=0.40 | **Yes — extremely** (p.adj=6.1×10⁻²⁶⁹) |
| 7 | Mouse GSE233866, 6-OHDA | SNc vs VTA, **within lesioned cells only** (1 combined network) | **Formal DME test** | blue, kME=0.32 | **Yes — extremely, and highly specific** (p.adj=4.5×10⁻¹²³, higher in SNc) |

Study #7 answers a different question than #3/#4 or #5/#6: not baseline regional structure, and
not disease-vs-healthy within one region, but **what's specifically more active in the vulnerable
region (SNc) than the resistant one (VTA) while neurodegeneration is actually happening** — the
comparison most directly relevant to a region-selective treatment target. 11 of 12 modules differ
significantly here (expected — SNc and VTA are fundamentally different populations), but
CACNA1D's own module is among the strongest signals in the whole table, and specifically more
active in SNc — arguably the single most direct result this project has produced for "why SNc
specifically."

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
match.

**Results — region-selective signal during neurodegeneration (the sharpest single finding)** —
study #7. This is the comparison most directly relevant to a region-selective treatment: within
lesioned tissue only, CACNA1D's module is one of the most significant differences in the entire
SNc-vs-VTA table (p.adj=4.5×10⁻¹²³) and is specifically *more active in SNc*, the region that's
actually dying. Lead with this alongside #6 as the strongest quantitative results in the whole
project — but note the caveat that 11 of 12 modules differ here, since SNc and VTA are different
populations to begin with, so frame CACNA1D's result as "notably strong even against that broad
background," not as a uniquely surprising isolated hit.

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
substitute for direct human validation of any proposed target.

## File map

- `human_GSE243639/pd_vs_control_separate_networks/{control_network,pd_network}/` — study #1
- `human_GSE243639/pd_vs_control_combined_dme/` — study #2
- `mouse_GSE233866/healthy_baseline/{snc_network,vta_network}/` — study #3
- `mouse_GSE233866/healthy_baseline/snc_vs_vta_preservation/` — study #4
- `mouse_GSE233866/lesion_model/{snc_lesioned_network,snc_intact_network}/` — study #5
- `mouse_GSE233866/lesion_model/lesioned_vs_intact_combined_dme_snc/` — study #6
- `mouse_GSE233866/lesion_model/vta_lesioned_network/` + `snc_vs_vta_combined_dme_lesioned/` — study #7
- `reminders/reminder.md` — running project log of decisions and pivots made along the way
