# hdWGCNA End-to-End System Check — Master Validation Prompt

> **How to use this file:** Paste this entire document into your AI coding assistant
> (Claude Code, Cursor, etc.) with your project repository open. It instructs the AI
> to perform a systematic, top-to-bottom audit of an hdWGCNA analysis pipeline — from
> environment and data integrity through every pipeline stage, the R microservice
> bridge, the Python backend, and the React frontend — and to report findings without
> silently "fixing" things in ways that could invalidate your science. Fill in the
> bracketed `[...]` placeholders with your project specifics before running.

---

## 0. Role and Objective (read first)

You are acting as a **senior computational biologist and research software engineer**
conducting a pre-submission validation audit of a Parkinson's disease single-nucleus
RNA-seq co-expression network analysis built on **hdWGCNA**. Your audience is a student
researcher who needs the pipeline to be scientifically correct, reproducible, and
defensible under peer review. The stakes are high: a silent bug in normalization,
grouping, or module assignment can produce a plausible-looking but wrong biological
conclusion. Your job is to find those bugs, not to make the code merely *run*.

**Prime directive:** Never trade scientific validity for the appearance of success.
A pipeline that errors loudly is safer than one that silently produces a corrupted
module table. When you find something wrong, report it; do not paper over it.

### How you must behave during this audit

1. **Verify, don't assume.** For every hdWGCNA/Seurat function call, confirm the
   argument names and defaults against the installed package version — not from memory.
   If you cannot verify a signature, say so explicitly and flag it as an open risk.
2. **State assumptions out loud.** Whenever you must assume something (a metadata
   column name, a file format, an expected cell count), write the assumption down in
   the report and mark what would happen if it's wrong.
3. **Identify data gaps.** If a check requires information you don't have (e.g., the
   number of donors, the sequencing platform), list it as a gap rather than inventing
   a value.
4. **Never fabricate results.** Do not invent gene names, module colors, kME values,
   p-values, enrichment terms, or citation details. If a value should come from
   running code, run the code or mark it as "requires execution."
5. **Distinguish severity.** Label every finding as `BLOCKER`, `WARNING`, or `NOTE`.
6. **Do not auto-edit scientific parameters.** You may fix syntax, imports, and
   obvious typos. You may NOT change soft power, `k`, `fraction`, `group.by`, network
   type, or any parameter that changes the biology without explicitly asking first.

Produce a single structured **Validation Report** at the end (template in Section 14).

---

## 1. Project Context (fill this in)

- **Research question:** Which genes are co-expressed with the L-type calcium channel
  genes `CACNA1D` (Cav1.3) and `CACNA1C` (Cav1.2) in substantia nigra pars compacta
  (SNc) dopaminergic neurons, and how do those co-expression modules relate to
  oxidative stress, mitochondrial/lysosomal dysfunction, and calcium-signaling
  pathways, comparing Parkinson's disease (PD) donors to healthy controls (HC).
- **Biological hypothesis under test:** The "calcium channel hypothesis" of selective
  SNc dopaminergic vulnerability (autonomous pacemaking via Cav1.3/Cav1.2 → chronic
  Ca²⁺ influx → mitochondrial oxidative stress).
- **Dataset(s):** `[e.g., Kamath et al. 2022 human midbrain snRNA-seq; Smajić et al.
  2022]`
- **Cell type of interest / `group_name`:** `[e.g., "DA" or "dopaminergic_neuron"]`
- **Condition column:** `[e.g., "disease" with levels "PD" / "HC"]`
- **Sample/donor column:** `[e.g., "donor_id" or "Sample"]`
- **hdWGCNA version target:** `[e.g., 0.4.12]`
- **Seurat version:** `[e.g., v5]`
- **Compute environment:** `[e.g., Docker r-service on Apple Silicon Mac, 16 GB RAM]`

---

## 2. Environment & Dependency Audit

Check the R and system environment before touching any data.

**2.1 R version and platform.** Confirm R ≥ 4.2 (hdWGCNA 0.4.x targets modern R;
the tutorial was compiled under R 4.4). Report `R.version.string` and platform.
On Apple Silicon, flag whether R is running natively (`aarch64`) or under Rosetta
(`x86_64`), because some Bioconductor binaries differ and WGCNA's compiled code is
sensitive to the BLAS in use.

**2.2 Package presence and versions.** Confirm each of these loads without error and
report its version: `Seurat`, `SeuratObject`, `hdWGCNA`, `WGCNA`, `igraph`,
`tidyverse`, `cowplot`, `patchwork`, `ggraph`, `UCell`, `GenomicRanges`,
`GeneOverlap`, `harmony`, `enrichR`. For each, note if the installed major version
could change behavior (e.g., Seurat v4 vs v5 assay/layer semantics).

**2.3 Seurat v4 vs v5 layer semantics.** This is a frequent silent failure. In Seurat
v5, expression data lives in *layers* (`counts`, `data`, `scale.data`) rather than
v4 *slots*. Confirm the code uses the correct accessor for the installed version and
that `SetDatExpr(..., layer = 'data')` (v5) vs `slot = 'data'` (v4) matches. Flag any
mismatch as a `BLOCKER` because it can silently pull the wrong matrix.

**2.4 WGCNA threading.** If `enableWGCNAThreads()` / `allowWGCNAThreads()` is used,
confirm it does not conflict with the container's CPU limits. Over-threading inside a
memory-constrained Docker container is a common cause of silent kills (OOM). Report
the thread count and available cores/RAM.

**2.5 Random seed.** Confirm a fixed seed (`set.seed(...)`) is set before any
stochastic step (metacell KNN, network construction). Without it, results are not
reproducible run-to-run. Mark missing seed as a `WARNING`.

**2.6 BLAS/LAPACK.** Report which BLAS is linked (`sessionInfo()`), since numerical
differences in correlation/eigendecomposition across BLAS implementations can subtly
shift module boundaries. Not a blocker, but note it for reproducibility.

---

## 3. Input Data Integrity Audit

Before any hdWGCNA function runs, validate the Seurat object itself.

**3.1 File loads.** Confirm the `.h5Seurat`/`.rds` path exists on the shared volume,
loads without error, and yields a Seurat object. If the loader is `SeuratDisk`,
confirm `Convert()`/`LoadH5Seurat()` completed and no assay was dropped.

**3.2 Object shape.** Report: number of cells, number of genes, assay names, default
assay, available reductions (`pca`, `harmony`, `umap`), and layers/slots present.
Confirm a dimensionality reduction suitable for metacell KNN exists (the tutorial
uses `harmony`).

**3.3 Metadata columns exist.** Confirm the exact column names referenced downstream
actually exist in `seurat_obj@meta.data`: the cell-type column, the condition column,
and the sample/donor column. Case-sensitive. A typo here (e.g., `cell_typea` vs
`cell_type`) is a `BLOCKER`. List all metadata column names verbatim so mismatches
are visible.

**3.4 Cell-type labels.** Print the unique values of the cell-type column. Confirm the
target `group_name` (dopaminergic neuron label) is present **exactly** as spelled in
the code. Report how many cells carry that label. If SNc DA neurons are <200 cells or
<a few per donor, flag that metacell construction may be unreliable for this cell type
(hdWGCNA documents poor behavior for extremely underrepresented cell types).

**3.5 Condition and donor structure.** Cross-tabulate condition × donor: how many PD
donors, how many HC donors, how many cells of the target type per donor per condition.
This matters enormously for the PD-vs-HC comparison downstream — a module "difference"
driven by a single donor is not a disease effect. Report the table and flag any
condition with only 1–2 donors as a `WARNING` (confounded, underpowered).

**3.6 Normalization state.** Confirm whether the `data` layer is log-normalized
(hdWGCNA expects normalized data for `SetDatExpr(layer='data')`). Check that counts
and normalized data are not accidentally swapped. Report summary statistics
(min/median/max) of a few housekeeping genes to sanity-check the scale.

**3.7 Target genes present.** Confirm `CACNA1D` and `CACNA1C` actually exist in the
gene set of the object (correct symbol, correct species — human). If the object uses
Ensembl IDs instead of symbols, flag that the gene-selection and module lookup must
use IDs. Report expression fraction of both genes in the target cell type — if either
is expressed in <5% of DA neurons, note that the `fraction` gene-selection threshold
may exclude your gene of interest entirely (a silent, project-ending bug).

---

## 4. Stage 1 — `SetupForWGCNA`

**4.1 Parameters.** Confirm `gene_select`, `fraction`, and `wgcna_name` are set and
sensible. If `gene_select = "fraction"` with `fraction = 0.05`, confirm that
`CACNA1D`/`CACNA1C` survive this filter in the target cell type (tie back to 3.7). If
they don't, the entire analysis cannot answer the research question — `BLOCKER`.

**4.2 Gene count after selection.** Report how many genes passed selection. Too few
(<1000) risks an unstable network; too many (>10,000) risks noise and long runtimes.
Note the number and whether it's in a reasonable band.

**4.3 Idempotency / no subsetting after setup.** hdWGCNA explicitly does **not**
support subsetting the Seurat object after `SetupForWGCNA`. Confirm the code does not
`subset()` the object at any point after this call. Flag any post-setup subsetting as
a `BLOCKER`.

---

## 5. Stage 2 — `MetacellsByGroups` + `NormalizeMetacells`

**5.1 `group.by` correctness.** Confirm `group.by` includes **both** the cell-type
column and the sample/donor column (tutorial: `c("cell_type", "Sample")`). This is the
single most important scientific check in metacell construction: if you do not group
by donor, you will aggregate cells across patients and manufacture fake co-expression.
Missing donor grouping is a `BLOCKER`.

**5.2 `reduction`.** Confirm the KNN reduction (e.g., `harmony`) exists and is
appropriate (batch-corrected space preferred so metacells aren't split by batch).

**5.3 `k` and `max_shared`.** Report both. Confirm `k` is in a defensible range
(commonly 20–75; smaller for smaller datasets). Confirm `max_shared` limits metacell
overlap so metacells aren't near-duplicates. Note the values but do NOT change them
without asking.

**5.4 `min_cells`.** Confirm groups smaller than `min_cells` are excluded, and that
this doesn't silently drop your PD or HC DA-neuron populations for some donors. Report
which groups were dropped, if any.

**5.5 Metacell yield.** After construction, report the number of metacells per
condition and per donor for the target cell type. If a condition ends up with very few
metacells, the downstream DME comparison is underpowered — `WARNING`.

**5.6 Normalization.** Confirm `NormalizeMetacells` ran and that downstream steps use
the normalized metacell matrix, not raw. Report the layer used.

---

## 6. Stage 3 — `SetDatExpr`

**6.1 Group selection.** Confirm `group_name` matches the DA-neuron label exactly and
`group.by` matches the column used in `MetacellsByGroups` (the tutorial warns these
must be the same column). Mismatch is a `BLOCKER`.

**6.2 Assay/layer.** Confirm `assay` and `layer`/`slot` pull normalized expression
(`layer='data'` in v5). Re-verify against 2.3.

**6.3 Multi-group intent.** If the design intends to build the network on DA neurons
only, confirm a single `group_name` is passed (not a vector), unless combining cell
types is deliberate. Note the choice.

**6.4 Expression matrix shape.** Report the dimensions of the stored expression matrix
(metacells × genes). Confirm it's transposed as hdWGCNA expects and contains no
all-zero or all-NA genes that would break correlation.

---

## 7. Stage 4 — `TestSoftPowers` + soft-power selection

**7.1 Network type consistency.** Confirm `networkType` (`signed` / `unsigned` /
`signed hybrid`) is set once and used consistently through `ConstructNetwork`. A
signed network in `TestSoftPowers` but unsigned in `ConstructNetwork` is a silent
inconsistency — `BLOCKER`.

**7.2 Power sweep ran.** Confirm `TestSoftPowers` produced a power table
(`GetPowerTable`). Report the table.

**7.3 Chosen power justified.** Confirm the selected soft power is the lowest power
achieving a scale-free topology fit (SFT.R.sq) ≥ ~0.8 (tutorial guidance). If the code
hard-codes a power, confirm it matches what the sweep supports for *this* dataset — do
not assume 9 (the tutorial's value) transfers. If no power reaches 0.8, flag it: the
data may be too sparse or the cell population too small — `WARNING`, and note that
`ConstructNetwork` will otherwise auto-pick a power.

**7.4 Reproducibility.** Confirm the power choice is recorded (logged or written to
the output directory) so the run is auditable.

---

## 8. Stage 5 — `ConstructNetwork`

**8.1 Soft power passed.** Confirm the justified power from Section 7 is actually
passed to `ConstructNetwork` (or that auto-selection is intentional and logged).

**8.2 TOM written.** Confirm `tom_name` is set and the TOM file is written to disk in
the output directory. Report the path and file size.

**8.3 Dendrogram sanity.** Confirm `PlotDendrogram` produces a figure with more than
one non-grey module. A network where nearly everything is grey (unassigned) indicates
a failed or degenerate network — `WARNING`/`BLOCKER` depending on severity.

**8.4 Grey module handling.** Confirm downstream code **excludes the grey module**
(unassigned genes) from interpretation, as the tutorial instructs. Including grey in
enrichment or DME is a scientific error — `BLOCKER` if present.

**8.5 Module count.** Report the number of modules. Very few (1–2) or very many (>50)
can both indicate poor parameterization; note for interpretation.

---

## 9. Stage 6 — `ModuleEigengenes` + `ModuleConnectivity`

**9.1 ScaleData prerequisite.** Confirm `ScaleData` (or the metacell scaling wrapper)
ran before `ModuleEigengenes` if Harmony correction is used, or the harmony step in
`ModuleEigengenes` will error. Report whether harmonized MEs (hMEs) are being computed
(`group.by.vars`).

**9.2 Harmonization column.** Confirm `group.by.vars` in `ModuleEigengenes` is the
sample/donor column — so technical between-donor variation is removed from the
eigengenes. This is important for a clean PD-vs-HC comparison. Note the column used.

**9.3 kME computed in the right group.** Confirm `ModuleConnectivity` is run with
`group.by`/`group_name` set to the DA-neuron population used for `ConstructNetwork`
(tutorial recommends computing kME in the same group). Mismatch weakens hub-gene
validity — `WARNING`.

**9.4 Module of interest.** Confirm the code locates which module contains `CACNA1D`
and `CACNA1C` (`GetModules()` then filter). Report the module color/name for each.
This is the pivot of the entire project — if this lookup is missing, the pipeline
does not actually answer the research question — `BLOCKER`.

**9.5 Hub genes.** Confirm `GetHubGenes` / `PlotKMEs` extracts top hub genes by kME.
Report whether the code captures hub genes for the CACNA1D/CACNA1C module specifically.

**9.6 Output persistence.** Confirm the final object is saved (`saveRDS`) and the
module assignment table is written to the output directory as CSV for the Python
backend to serve.

---

## 10. Downstream Scientific Validity (the part that determines whether the paper holds)

**10.1 Co-localization check.** Do `CACNA1D` and `CACNA1C` fall in the same module or
closely correlated modules? Report the answer and treat either outcome as a finding,
not a failure. Do not tune parameters to force them together — that would be
circular reasoning. Flag any code that appears to do so as a `BLOCKER`.

**10.2 Enrichment sanity.** If enrichment (`RunEnrichr` or `enrichR`) is run on the
CACNA1D module, confirm the gene background is correct (the WGCNA gene set, not the
whole genome) and that the databases queried are appropriate (GO Biological Process,
KEGG, Reactome). Confirm the code does **not** cherry-pick only calcium/mitochondrial
terms while hiding others — report the full top-N enrichment, favorable or not.

**10.3 PD vs HC comparison design.** If `FindDMEs` / differential ME analysis compares
PD vs HC, confirm the comparison is at the correct level and controls for donor. A
difference driven by cell count imbalance or a single donor is not a disease effect.
Confirm the test and grouping. Flag pseudoreplication (treating cells as independent
when they're nested within a few donors) as a `WARNING` — this is one of the most
common and serious errors in single-cell differential analysis.

**10.4 Multiple testing.** If many modules or genes are tested for PD-vs-HC
differences, confirm p-values are corrected (FDR/Benjamini-Hochberg). Uncorrected
p-values across dozens of modules will produce false positives — `WARNING`.

**10.5 Direction and effect size.** Confirm reported differences include direction
(up/down in PD) and an effect-size measure, not just significance. Note if missing.

**10.6 Negative-control awareness.** Suggest (do not fabricate) at least one sanity
check the researcher could run: e.g., confirm a housekeeping module shows no PD-vs-HC
difference, or that a known PD gene (`SNCA`, `LRRK2`, `PINK1`, `PRKN`) behaves as
expected. Flag if no such control exists in the analysis.

---

## 11. R Microservice (plumber bridge) Audit

**11.1 Endpoint ↔ function fidelity.** Confirm each plumber endpoint calls the real
hdWGCNA/Seurat function with correct argument names (re-verify against installed
package). No reimplementation, no silent parameter substitution.

**11.2 Parameter passthrough.** Confirm the values entered in the UI (h5Seurat path,
output dir, cell-type column, `group.by[0]`, `group.by[1]`, group name, `wgcna_name`,
network type, `k`, `max_shared`) are passed through unchanged to the R functions.
Trace one full request end-to-end and confirm no field is dropped, reordered, or
defaulted silently. The UI's `GROUP BY [0]`/`GROUP BY [1]` must map to the
`group.by = c(col0, col1)` vector in `MetacellsByGroups` — confirm order and names.

**11.3 Path safety.** Confirm the service validates that the h5Seurat path exists and
is inside the shared volume before running, returning a clear error otherwise (not a
generic 500). Confirm the output directory is created if absent.

**11.4 Async/job handling.** If jobs run via `future`/async, confirm job status is
tracked, errors propagate to the status endpoint (not swallowed), and a failed R step
sets status to `error` with a readable message rather than hanging on `running`.

**11.5 Result serialization.** Confirm the module table is serialized to JSON/CSV
faithfully (numeric precision for kME preserved) and plots are written as PNG to the
output directory and served correctly.

**11.6 Resource limits.** Confirm the container has enough RAM for the dataset;
network construction on large gene sets is memory-heavy. Flag if the Mac's Docker VM
memory limit is likely to OOM-kill the job.

---

## 12. Python Backend Audit

**12.1 Contract match.** Confirm the FastAPI request/response schema matches what the
R service expects and what the React UI sends. Field names, types, and required-ness
should line up across all three layers. Report any drift.

**12.2 File handling.** Confirm uploaded files land on the shared volume that the R
container also mounts, and that the path handed to R is the path *inside the R
container*, not the Python container (a classic cross-container path bug). `BLOCKER`
if the two services don't agree on the path.

**12.3 Polling/status.** Confirm the backend correctly polls the R job status and does
not report `SUCCESS` before results exist. Confirm error states surface to the user.

**12.4 Timeouts.** Confirm HTTP timeouts are long enough for real network construction
(minutes to tens of minutes), or that the design is fully async so the request doesn't
time out mid-analysis.

**12.5 No hidden reprocessing.** Confirm the Python layer does not re-normalize,
re-filter, or otherwise mutate the expression data — it should only orchestrate, not
do science. Any silent data transformation in Python is a `BLOCKER`.

---

## 13. React Frontend Audit

**13.1 Field validation.** Confirm the form validates required fields (non-empty
h5Seurat path, valid cell-type column) before submitting, so a stray value like a
single character `O` in the path can't launch a doomed job. Report current validation.

**13.2 Value fidelity.** Confirm the exact strings typed by the user reach the backend
unchanged (no trimming that would break a legitimate label, no autocomplete
corruption). The `cell_typea` typo visible in the UI is exactly the class of error
front-end validation should catch — recommend validating the cell-type column against
the actual metadata columns returned from the object.

**13.3 State/status display.** Confirm the UI reflects real job status (queued /
running / done / error) from the backend and doesn't show stale or optimistic states.

**13.4 Result rendering.** Confirm module tables and plots render from the actual
result endpoints and are tied to the correct `job_id`.

**13.5 No browser storage misuse.** Confirm no reliance on unsupported browser storage
in a way that would break the app in its host environment.

---

## 14. Required Output — Validation Report Template

Produce the report in exactly this structure:

```
# hdWGCNA Pipeline Validation Report
Date: <date>    Reviewer: AI system-check    Repo commit: <hash if available>

## Executive summary
- Overall status: PASS / PASS-WITH-WARNINGS / FAIL
- Count of findings: X BLOCKER, Y WARNING, Z NOTE
- One-paragraph plain-language verdict on whether results are trustworthy.

## Verified environment
- R version, platform (native vs Rosetta), key package versions, seed set (y/n), BLAS.

## Findings (ordered: all BLOCKERs first, then WARNINGs, then NOTEs)
For each finding:
- ID: <e.g., B1, W3, N2>
- Severity: BLOCKER | WARNING | NOTE
- Location: <file:line or pipeline stage>
- What's wrong:
- Why it matters (scientific or engineering consequence):
- Evidence (what you checked / observed):
- Recommended fix (do NOT auto-apply if it changes biology):

## Assumptions I had to make
- List every assumption, and what breaks if it's false.

## Data gaps (information I needed but did not have)
- List each gap and which check it blocked.

## Scientific validity summary
- CACNA1D module: <color/name or "not determined">
- CACNA1C module: <color/name or "not determined">
- Co-localized? <yes/no/uncertain>
- PD-vs-HC comparison sound? <assessment, incl. donor-level replication>
- Multiple-testing correction present? <y/n>

## What I did NOT change (and why)
- List parameters/biology I intentionally left for the researcher to decide.

## Suggested next validation steps for the researcher
- Concrete, non-fabricated checks to run manually.
```

At the very end of your report, do a final self-check and answer these four questions
explicitly: (1) Did I verify function signatures against the installed package rather
than memory? (2) Did I avoid fabricating any gene, module, statistic, or citation?
(3) Did I flag rather than silently fix anything that changes the biology? (4) Did I
clearly separate what I verified from what I assumed?

---

## 15. Guardrails (do not violate)

- Do **not** reimplement or approximate any hdWGCNA/WGCNA algorithm.
- Do **not** change scientific parameters (soft power, `k`, `fraction`, `group.by`,
  network type, `max_shared`) to make something pass — flag and ask.
- Do **not** fabricate gene names, module assignments, kME values, p-values,
  enrichment terms, dataset statistics, or citations.
- Do **not** tune parameters to force `CACNA1D`/`CACNA1C` into the same module.
- Do **not** report `SUCCESS` for any layer whose outputs you could not actually
  verify — say "requires execution" instead.
- Do **not** hide unfavorable enrichment or differential results.
- When uncertain, downgrade confidence and say so. A hedged, honest report is the
  deliverable — not a green checkmark.

---

## 16. Scientific Background, Credits & Citations

### Project scientific framing
This pipeline investigates the **calcium channel hypothesis** of selective dopaminergic
neuron vulnerability in Parkinson's disease, associated with **D. James Surmeier and
colleagues**, which proposes that substantia nigra pars compacta dopaminergic neurons
rely on Cav1.3 (`CACNA1D`) and Cav1.2 (`CACNA1C`) L-type calcium channels for
autonomous pacemaking, producing chronic Ca²⁺ influx, elevated mitochondrial oxidative
stress, and selective vulnerability. The dihydropyridine blocker **isradipine** was
tested clinically on this basis and did not meet its primary endpoint in the Phase 3
**STEADY-PD III** trial, motivating co-expression-network approaches that examine genes
co-regulated with these channels rather than the channels in isolation.

> Note: the citations below are provided as starting references to verify and format in
> your reference manager. Confirm every detail (authors, year, volume, DOI) against the
> primary source before including in a manuscript — do not treat this list as
> authoritative without checking.

### Software to cite

**hdWGCNA (primary method)**
Morabito S, Reese F, Rahimzadeh N, Miyoshi E, Swarup V. "hdWGCNA identifies
co-expression networks in high-dimensional transcriptomics data." *Cell Reports
Methods*, 2023. Package: https://github.com/smorabit/hdWGCNA ; documentation:
https://smorabit.github.io/hdWGCNA/. Developed by Sam Morabito, Swarup Lab, UC Irvine.

**WGCNA (underlying algorithm)**
Langfelder P, Horvath S. "WGCNA: an R package for weighted correlation network
analysis." *BMC Bioinformatics*, 2008;9:559.
Zhang B, Horvath S. "A general framework for weighted gene co-expression network
analysis." *Statistical Applications in Genetics and Molecular Biology*, 2005;4:Article17.

**Seurat (single-cell framework)**
Hao Y, et al. "Dictionary learning for integrative, multimodal, and scalable
single-cell analysis." *Nature Biotechnology*, 2024 (Seurat v5). Plus earlier Seurat
papers (Butler et al. 2018; Stuart et al. 2019; Hao et al. 2021) as appropriate to the
version used.

**Harmony (batch integration)**
Korsunsky I, et al. "Fast, sensitive and accurate integration of single-cell data with
Harmony." *Nature Methods*, 2019;16:1289–1296.

**UCell (module gene scoring)**
Andreatta M, Carmona SJ. "UCell: Robust and scalable single-cell gene signature
scoring." *Computational and Structural Biotechnology Journal*, 2021.

**Enrichr / enrichR (enrichment)**
Chen EY, et al. "Enrichr: interactive and collaborative HTML5 gene list enrichment
analysis tool." *BMC Bioinformatics*, 2013. Kuleshov MV, et al. *Nucleic Acids
Research*, 2016.

**Supporting R ecosystem**
R Core Team (R statistical computing environment); `tidyverse` (Wickham et al., *JOSS*
2019); `igraph` (Csárdi & Nepusz, 2006); `ggraph`/`patchwork`/`cowplot` (Pedersen;
Wilke).

### Datasets to cite (confirm the one(s) you actually use)

**Kamath et al. 2022** — "Single-cell genomic profiling of human dopamine neurons
identifies a population that selectively degenerates in Parkinson's disease."
*Nature Neuroscience*, 2022;25:588–595.

**Smajić et al. 2022** — "Single-cell sequencing of human midbrain reveals glial
activation and a Parkinson-specific neuronal state." *Brain*, 2022;145(3):964–978.

### Scientific hypothesis references (verify before citing)

Representative Surmeier-group and related work on the calcium/pacemaking/vulnerability
hypothesis (confirm exact articles and years against the primary literature):
Chan CS, et al. *Nature*, 2007 (rejuvenation of SNc neurons via Cav1.3). Guzman JN,
et al. *Nature*, 2010 (oxidant stress from pacemaking). Surmeier DJ, et al. — reviews
on calcium and mitochondrial stress in PD vulnerability.

**STEADY-PD III trial**
Parkinson Study Group STEADY-PD III Investigators. "Isradipine versus placebo in early
Parkinson disease: a randomized trial." *Annals of Internal Medicine*, 2020. Confirm
authorship, volume, and pages before citing.

### Credits
- **hdWGCNA method & software:** Sam Morabito and the Swarup Lab, University of
  California, Irvine.
- **WGCNA framework:** Peter Langfelder and Steve Horvath, UCLA.
- **Calcium channel hypothesis:** D. James Surmeier and colleagues.
- **Pipeline assembly, microservice bridge, and application code:** `[your name /
  lab / affiliation]`.
- **Validation prompt / system-check design:** prepared as a research support artifact
  for this project.

---

---

## Appendix A — Per-Stage Diagnostic Snippets (verify, then run)

These are reference diagnostics the AI can adapt to the actual object and package
version. **Do not run blindly** — confirm each function/argument exists in the
installed version first, and adjust column names to match the real metadata. They are
provided so the audit produces concrete evidence, not vibes.

**A.1 Environment fingerprint**
```r
R.version.string
sessionInfo()                      # BLAS/LAPACK, attached package versions
packageVersion("hdWGCNA"); packageVersion("Seurat"); packageVersion("WGCNA")
Sys.getenv("OMP_NUM_THREADS")      # threading sanity in container
```

**A.2 Object integrity**
```r
seurat_obj
dim(seurat_obj)                    # genes x cells
Assays(seurat_obj); DefaultAssay(seurat_obj)
Reductions(seurat_obj)             # expect pca / harmony / umap
colnames(seurat_obj@meta.data)     # confirm exact metadata column names
```

**A.3 Metadata / grouping sanity**
```r
table(seurat_obj$cell_type)                       # unique cell-type labels + counts
table(seurat_obj$disease)                          # condition levels (PD/HC)
table(seurat_obj$cell_type, seurat_obj$disease)    # target type present in both?
with(seurat_obj@meta.data,
     table(donor_id, disease))                     # donors per condition (power!)
```

**A.4 Target gene presence & expression**
```r
c("CACNA1D","CACNA1C") %in% rownames(seurat_obj)   # both TRUE?
# fraction of DA neurons expressing each gene (tie to the `fraction` threshold):
da <- subset(seurat_obj, cell_type == "DA")        # inspection only — do NOT
expr <- GetAssayData(da, layer = "data")           #   subset the WGCNA object itself
rowMeans(expr[c("CACNA1D","CACNA1C"), ] > 0)
```

**A.5 Soft power evidence**
```r
GetPowerTable(seurat_obj)                          # SFT.R.sq column drives the choice
# lowest Power with SFT.R.sq >= ~0.8 is the defensible pick
```

**A.6 Module of interest lookup (the project pivot)**
```r
mods <- GetModules(seurat_obj)
mods[mods$gene_name %in% c("CACNA1D","CACNA1C"),
     c("gene_name","module")]                      # which module(s)?
table(mods$module)                                  # module sizes; watch grey
```

**A.7 Hub genes for the CACNA1D module**
```r
hubs <- GetHubGenes(seurat_obj, n_hubs = 25)
target_mod <- mods$module[mods$gene_name == "CACNA1D"]
hubs[hubs$module == target_mod, ]
```

---

## Appendix B — Healthy Result vs. Red Flag (quick reference)

| Stage | Looks healthy | Red flag (investigate before trusting results) |
|---|---|---|
| Environment | Native aarch64 R, seed set, versions logged | Rosetta mismatch, no seed, unknown versions |
| Data load | Object loads, reductions present | Missing harmony/pca, dropped assay |
| Metadata | Target labels + PD/HC both present, ≥3 donors each | Typo'd column, one condition has 1 donor |
| Target genes | CACNA1D/CACNA1C present, expressed >5% of DA | Absent, or below the `fraction` cutoff |
| Setup | 1,000–10,000 genes selected | <500 or >15,000 genes, or post-setup subsetting |
| Metacells | group.by includes donor; sensible k; both conditions yield metacells | No donor grouping; a condition drops out |
| SetDatExpr | Same column as metacells; normalized layer | Column mismatch; raw counts pulled |
| Soft power | Lowest power with SFT ≥0.8, logged | Hard-coded 9 with no check; no power hits 0.8 |
| Network | >1 non-grey module; TOM written | Nearly all grey; degenerate network |
| Eigengenes | hMEs harmonized by donor; kME in DA group | Harmonized by wrong column; kME in wrong group |
| Module lookup | CACNA1D module identified | Lookup step absent |
| PD vs HC | Donor-level replication; FDR-corrected; effect sizes | Pseudoreplication; uncorrected p; no direction |

---

## Appendix C — The Five Most Likely Silent Bugs in This Specific Project

Ranked by how often they cause plausible-but-wrong results, and how hard they are to
notice after the fact:

1. **Metacells not grouped by donor.** Manufactures cross-patient co-expression. The
   network will look great and mean nothing. Check `group.by` first, always.
2. **CACNA1D filtered out by the `fraction` threshold.** L-type channels can be
   sparsely detected in snRNA-seq; if `CACNA1D` is expressed in <5% of DA neurons it
   silently vanishes at `SetupForWGCNA`, and the whole project can't answer its
   question. Verify presence *after* gene selection, not just in the raw object.
3. **Seurat v4/v5 layer mismatch.** Pulls the wrong matrix (raw vs normalized) into
   `SetDatExpr`, corrupting every correlation downstream with no error thrown.
4. **Pseudoreplication in PD-vs-HC.** Treating thousands of cells as independent when
   they come from a handful of donors inflates significance massively. A "highly
   significant" module difference may be one unusual donor.
5. **Cross-container path mismatch.** The Python container's file path is handed to the
   R container, which can't see it; job fails confusingly, or worse, silently reads a
   stale/empty file if a path happens to resolve.

---

*End of master validation prompt. Fill in bracketed fields, paste into your AI coding
assistant with the repository open, and require it to return the Section 14 report.
Confirm all citations against primary sources before using them in any manuscript.*
