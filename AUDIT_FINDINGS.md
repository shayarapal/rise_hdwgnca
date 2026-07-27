# hdWGCNA Pipeline — Validation Audit: Running Issues Log

Live findings from the `HDWGCNA_SYSTEM_CHECK_PROMPT.md` audit, accumulated step by step.
Rolls up into the Section 14 Validation Report when the audit is complete.

- **Target object:** `DATA_UNZIPPED/seurat_GSE243639_SNc_DA.rds` (GSE243639, SNc DA neurons)
- **Verification env:** native conda `hdWGCNA` (R 4.5.3, hdWGCNA 0.4.12, Seurat 5.5.0) on this Mac
- **Production env (PC):** Windows + Docker Compose — *not yet fingerprinted (build in progress)*
- **Research pivot:** co-expression partners of `CACNA1D` (Cav1.3) + `CACNA1C` (Cav1.2), PD vs control

**Severity:** `BLOCKER` (invalidates results) · `WARNING` (threatens validity/reproducibility) · `NOTE` (record only)
**Status:** OPEN · VERIFYING · RESOLVED · ACCEPTED (user-acknowledged, no change) · DEFERRED (checked in a later step)

---

## Master issue table

| ID | Sev | Step/§ | Title | Status |
|----|-----|--------|-------|--------|
| **B1** | **BLOCKER** | §11 / Docker | **r-service Docker image silently ships WITHOUT hdWGCNA/WGCNA** — build exits 0 but pipeline cannot run. WGCNA's Bioconductor deps (GO.db, impute, preprocessCore, AnnotationDbi) missing from Dockerfile; same latent bug in `install.R` | **RESOLVED** — rebuilt image verified: hdWGCNA 0.4.12 + WGCNA 1.73 + all Bioc deps present; build tripwire added |
| **B2** | **BLOCKER (scientific, parameter-conditional)** | 7 / §8 | At the auto-recommended power **4**: network is **77.5% grey** and **both `CACNA1D` and `CACNA1C` are unassigned (grey)** → the research question is **unanswerable as currently parameterized**. Not a code bug; needs soft-power re-evaluation | OPEN — methodological (do NOT tune to force targets together) |
| **B3** | **BLOCKER** | §11 / Docker | r-service Dockerfile `apt-get` missing **libsodium-dev** (→ `sodium` → **plumber** fails) and **libuv1-dev** (→ `fs`→sass→bslib→ **tidyverse** fails). Silently absent; container crash-loops at `library(plumber)`. Found by first actual container run | **RESOLVED** — rebuilt (all 14 runtime pkgs load); full Docker stack smoke-tested end-to-end |
| W1 | WARNING | 1 / §2 | Env drift: Docker is R **4.4.3** (native 4.5.3); moot until B1 fixed since the image can't run at all | SUPERSEDED by B1 |
| W2 | WARNING | 1 / §2 | Thread oversubscription: `OMP_NUM_THREADS` unset under WGCNA×future workers | ACCEPTED (user confident on PC hardware) |
| W3 | WARNING | 2 / §3 | `CACNA1C` survives `fraction=0.05` by a thin margin — drops if fraction > ~0.062 | OPEN |
| W4 | WARNING | 2 / §3 | PD arm underpowered: 4 PD vs 7 control donors clear `min_cells=100` | OPEN (inherent to disease; document) |
| W5 | WARNING | 4 / §5 | Metacell imbalance + donor concentration: control 3,888 vs PD 1,594; 3 control donors capped at `target_metacells=1000` → 3,000/3,888 control metacells from 3 donors | OPEN |
| N1 | NOTE | 1 / §2 | R 4.5.3 exceeds tutorial/Docker R 4.4 | OPEN |
| N2 | NOTE | 2 / §3 | `scale.data` layer dropped (memory fix) — it DID trip `ModuleEigengenes` (was the root cause of B4) | RESOLVED via B4 fix (re-scales in `build_network`) |
| N3 | NOTE | 3 / §4 | Selected gene set (7,465) is on the high end of the sane band | OPEN |
| N4 | NOTE | 4 / §5 | Documented tutorial deviations: `reduction="pca"` (not harmony), `max_shared=15` (tutorial 10); `target_metacells` not request-tunable | OPEN (deliberate per CHANGELOG) |
| N5 | NOTE | 6 / §7 | Soft power computed on POOLED metacells; for `/module-preservation` it must be recomputed on the reference condition (code supports this via `condition_col` on `/test-soft-powers`) | OPEN (correct usage required) |
| **B4** | **BLOCKER → RESOLVED** | 8 / §9 | Harmonized `ModuleEigengenes` failed `non-conformable arguments` because the DA object dropped `scale.data` (N2) but kept its `ScaleData` @command, so hdWGCNA's guard passed and harmonization ran on an empty matrix. **Fixed**: `build_network` now runs `ScaleData(features=GetWGCNAGenes())` before `ModuleEigengenes`. `/analyze` verified **done** end-to-end (first full run this session): `modules.csv` = 7,465 rows w/ kME cols + kME plot | **RESOLVED** |
| **W7** | **WARNING (RAM) → mitigated** | 8 / diag | Full `/analyze` peaks **~9.4 GB RSS** at `fraction=0.05`/7,465 genes (not the CHANGELOG's "0.66 GB" — that was soft-power only, no TOM). User upgraded to **15 GB usable** → single `/analyze` fits. **Caveat: `entrypoint.R` runs 2 future workers — 2 concurrent `/analyze` jobs ≈ 19 GB > 15 GB.** Keep to one heavy job at a time, or lower `workers` to 1 | OPEN (single-job OK; document concurrency limit) |
| N6 | NOTE | 8 / §9 | `write_module_outputs` never `saveRDS`s the network object (only `modules.csv`); `PlotKMEs` used but no `GetHubGenes` hub-gene table exported | OPEN |
| W8 | WARNING | 9 / §10 | Research question's downstream is **not implemented**: no enrichment (`enrichR` attached but never called), no differential ME test (`FindDMEs`), no trait correlation. PD-vs-control has ONLY module preservation. Pathway enrichment of the CACNA module — a core aim — cannot be produced by the pipeline | OPEN (scope; "planned" in docs/SKILLS.md) |
| N7 | NOTE | 9 / §10 | No built-in positive/negative control. Recommend a manual sanity check: DA-identity genes TH/SLC6A3/SLC18A2 (80–95% detected) should co-cluster; if they scatter, distrust the network | OPEN (recommendation) |

## Verified-GOOD (passed checks — no action)

- **§2** Seurat v5 layer semantics correct: `SetDatExpr(layer="data")` at plumber.R:272/560/568 matches installed Seurat 5.5.0.
- **§2** Reproducibility seeds present: `set.seed(42)` before both `MetacellsByGroups` (plumber.R:254,542); `ConstructNetwork` (randomSeed=12345) and `ModulePreservation` (seed=12345) seeded by default.
- **§2** Native run is native x86_64 (no Rosetta).
- **§3** Both target genes present as **symbols** and **confirmed in the 7,465-gene selected set** via real `SelectNetworkGenes`.
- **§3** Normalization correct: `data` log-normalized (max 7.97), distinct from integer `counts` (max 1154); not swapped.
- **§3** `harmony` reduction present; all referenced metadata columns exist (no typos).
- **§4** No `subset()` after `SetupForWGCNA` in any of the 4 job entrypoints — all condition subsetting precedes setup (clears the prompt's flagged BLOCKER risk).
- **§5** Metacells are **donor-pure** (verified): `group.by=c("cell_type","sample_id")` → every one of 5,482 metacells maps to a single `sample_id`, 11 distinct donors, zero cross-patient mixing. **Clears the prompt's #1 silent-bug risk.**
- **§5** `NormalizeMetacells` ran; metacell `data` layer is log-normalized (non-integer, max 5.78).
- **§6** `datExpr` clean: 5,482 metacells × 7,465 genes, correct WGCNA orientation (metacells in rows), and **0 all-NA / 0 any-NA / 0 all-zero / 0 zero-variance genes** — nothing will NA-out `bicor`.
- **§6** Metacell aggregation *rescues* the sparse target genes: `CACNA1D` 7.98%→**78%** of metacells, `CACNA1C` 6.18%→**71%**, both with real variance — robust correlation signal once included (context for W3).
- **§11 / Docker** Full stack smoke-tested end-to-end after B1+B3 fixes: all 3 containers run, `POST /gene-selection` via `:8200`→`:8100` returned correct results (7,465 genes @ 0.05 — **identical to native**, cross-env reproducible), bind mount round-trips both directions (host `.rds` in, artifacts out to `local_data/results/`), frontend `:4000` serves + proxies new routes. Windows+Docker deployment validated.

## Data gaps

- **D1:** PC core count — needed to finalize W2 thread settings. *(User: PC can handle it; treating W2 as accepted.)*
- **D2:** Docker r-service image's actual R / BLAS / resolved hdWGCNA version — closes W1. *Pending build.*
- **D3:** Whether target-gene selection differs on the full multi-cell-type object (pipeline uses the DA object, so not blocking).

---

## Detailed findings

### B1 (BLOCKER) — Docker r-service image is silently missing hdWGCNA — the PC deployment cannot run

*Surfaced by building the images (the portability step), confirmed by running R inside the built image.*

**Symptom:** `docker compose build r-service` exits **0**, but the image has neither `hdWGCNA` nor `WGCNA`. Every pipeline endpoint would fail at `library(hdWGCNA)` inside the future worker — the whole service is non-functional on the PC.

**Confirmed in-image:** `hdWGCNA FALSE, WGCNA FALSE, GO.db FALSE, impute FALSE, preprocessCore FALSE, AnnotationDbi FALSE` (UCell/GenomicRanges/Seurat/SeuratDisk TRUE). R version in image: **4.4.3**.

**Root cause:** `r-service/Dockerfile` lists `WGCNA` in the CRAN `install.packages()` call, but WGCNA's four **Bioconductor** dependencies — `GO.db, impute, preprocessCore, AnnotationDbi` — are never installed (the `BiocManager::install(...)` step only lists `BiocGenerics, GenomicRanges, GeneOverlap, UCell`). So WGCNA can't install from CRAN alone → hdWGCNA's `install_github` fails with `dependency 'WGCNA' is not available` → **but `remotes` only issues a Warning, so `Rscript` exits 0 and the Docker RUN "succeeds."** Classic silent failure.

**Why it went unnoticed:** every prior test this session ran against the **native conda env**, which already had WGCNA + all Bioconductor deps from conda-forge/bioconda. The Docker build path was never exercised until now — exactly why the plan built the image before shipping to the PC.

**Same bug in `r-service/install.R`** (native non-Docker install): it sets CRAN-only repos and its `BiocManager::install()` omits the same four deps, so a native install on the PC (without conda) would fail identically.

**Fix (ready — infrastructure, not science):**
1. In `r-service/Dockerfile`, add `GO.db, impute, preprocessCore, AnnotationDbi` to the `BiocManager::install(...)` list, and ensure that Bioc step runs **before** WGCNA is installed (either move WGCNA out of the first CRAN block to after the Bioc step, or install these deps first).
2. Mirror the same additions in `r-service/install.R`.
3. Add a build-time tripwire so this can never silently pass again: a final `RUN Rscript -e 'library(hdWGCNA); library(WGCNA)'` — it exits non-zero if either is missing, failing the build loudly.
4. Rebuild and re-fingerprint (also closes W1's R-version question).

### B4 (BLOCKER, confirmed both environments) — /analyze fails at harmonized ModuleEigengenes

The pipeline's `write_module_outputs` path runs `ModuleEigengenes(group.by.vars=p$group_by[[2]])` (= `"sample_id"`, plumber.R:368-370). This **errors**:

- Isolation (native, harmony 2.0.5): `ModuleEigengenes()` no-harmony → **OK**; `ModuleEigengenes(group.by.vars="sample_id")` → **ERROR `non-conformable arguments`** (+ warnings `Layer 'scale.data' is empty`, `Number of dimensions changing from 30 to 50`).
- End-to-end `/analyze` through **Docker (harmony 1.2.3)**: ran ~530s, then **`failed: non-conformable arguments`** — the same error. So it is **not** a harmony-version incompatibility.

**Never caught before because `/analyze` had never been run end-to-end** (only `/test-soft-powers`, which stops before eigengenes, and `/gene-selection`). This audit's Docker `/analyze` run is the first full pipeline execution on real data — and it fails.

**Most likely root cause (entangled with B2):** at power 4 the network has only ~5 non-grey modules, so the module-eigengene matrix is ~5-dimensional; harmony's internal PCA (the `30→50` dims message) can't operate on so few dimensions → non-conformable matrix multiply. If a higher, better-fitting power (B2 diagnosis) yields more modules, harmony may succeed — which would resolve **both** B2 and B4. The alternative resolutions (run without harmonization = lose donor batch correction; or a code workaround inside hdWGCNA's harmony call = forbidden reimplementation) are worse.

**Decision needed:** the B2 power diagnostic (test powers 6/9) is now double-motivated — it diagnoses whether B4's blocker is caused by low module count. Reported honestly either way; not tuning to force the targets together.

### B3 (BLOCKER) — Docker image crash-loops: plumber + tidyverse silently missing (system libs)

Found by the first *actual* container run (`docker compose up r-service`): the container restart-loops with `Error in library(plumber): there is no package called 'plumber'`. Full package fingerprint of the image: `plumber MISSING, tidyverse MISSING` (everything else, incl. the B1-fixed hdWGCNA/WGCNA, present).

**Root cause (build log):** the `apt-get` list omits two system libraries:
- `libsodium-dev` → the R `sodium` package fails to compile (`sodium.h: No such file`) → **`plumber` fails** (it depends on `sodium`).
- `libuv1-dev` → the R `fs` package fails to configure (`libuv was not found`) → cascades `fs`→`sass`→`bslib`→`rmarkdown` → **`tidyverse` fails**.

`install.packages()` skips a failed package with only a warning and exits 0, so the build "succeeded" with a broken image — the exact same silent-failure class as B1, one layer down. My B1 tripwire only loaded hdWGCNA/WGCNA, so it didn't catch the missing plumber stack. Both were latent in the *original* image too; only surfaced now because this is the first time the container was actually started (all prior testing ran natively via conda, which supplies these libs).

**Fix applied:**
1. Added `libsodium-dev` + `libuv1-dev` to the Dockerfile `apt-get`.
2. **Hardened the tripwire** to load all 14 runtime packages (entrypoint.R + plumber.R top-level + `pipeline_libs()`), so any future missing runtime dep fails the build loudly instead of shipping.
3. Full rebuild required (apt layer change invalidates the R-package cache) — running in background.

**Note:** `install.R` (native path) doesn't `apt-get`; a native non-Docker install relies on the OS/conda to provide libsodium+libuv (conda did, which is why native worked all session). Docker is the user's chosen path, so the fix lives in the Dockerfile.

### B2 (BLOCKER, scientific / parameter-conditional) — target genes land in grey at the recommended power

Built the real network at `soft_power=4` (the §7 auto-recommendation) on the DA object.

- **8.2 TOM written — PASS.** `audit_TOM.rda`, 202.9 MB, in the out_dir.
- **8.5 module count — OK on count.** 5 non-grey modules (turquoise 685, blue 478, brown 292, yellow 170, green 53) + grey. Not 1–2, not >50.
- **8.3/8.4 GREY FRACTION — FAIL.** **5,787 / 7,465 genes (77.5%) are grey.** Dendrogram merges at height 0.84–0.94 (weak modular structure). The audit calls a mostly-grey network degenerate.
- **Project pivot — FAIL.** `GetModules` puts **both `CACNA1D` and `CACNA1C` in grey** (unassigned). Since grey must be excluded from interpretation (§8.4), there is **no module to report co-expression partners from** — the entire research question cannot be answered from this network as built.

**Diagnosis (not yet a fix):** power 4 is the *smallest* power with SFT.R.sq ≥ 0.8, but its **mean connectivity is 495** — very high, i.e. the network is still dense/noisy there. Scale-free fit keeps climbing well past 4 (R² 0.837→0.957 at power 6→0.988 at power 9, mean.k 495→132→19). The "smallest ≥ 0.8" heuristic (which hdWGCNA's own `ConstructNetwork(soft_power=NULL)` also uses, and which the `/soft-powers` endpoint recommends) can land on a poorly-modularized power. This is a **methodological** issue, not a code defect — the recommendation faithfully implements the documented rule.

**What must NOT happen:** tuning the power (or `fraction`/`minModuleSize`) specifically to push `CACNA1D`/`CACNA1C` into a module would be circular reasoning (prompt guardrail). Any re-evaluation must report the outcome honestly, favorable or not.

**Power diagnostic RESULT (powers 4/6/9, memory-bounded run):**

| Power | Modules | Grey % | CACNA1D | CACNA1C | Harmonized ME |
|---|---|---|---|---|---|
| 4 | 5 | 77.5% | grey | grey | FAIL |
| 6 | 6 | 70.1% | grey | blue | FAIL |
| 9 | 6 | 75.4% | turquoise | turquoise | FAIL |

**Conclusions:**
- **B2 is inherent, not power-fixable:** grey stays 70–77% at every power — this sparse snRNA-seq DA population does not form tight co-expression modules.
- **Target co-localization is power-dependent and CANNOT be claimed:** apart at power 6 (grey vs blue), together at power 9 (both turquoise = the largest module). The "same module" answer exists only at power 9 — reporting it would be power-cherry-picking (guardrail violation). Honest negative result: the data does not robustly support Cav1.2/Cav1.3 co-localization. Both landing in the *largest* module at power 9 also smells like being swept into the dominant module, not specific co-expression.
- **B4 is confirmed code-level** (fails at all powers — see B4 row).

**Connects to:** N5 (power computed on pooled, control-dominated metacells — W5), and the general sparsity of the target genes (W3: 6–8% of single cells).

### Step 9 — §10 Downstream Scientific Validity

**10.1 Co-localization — see the power diagnostic (B2).** Power-unstable; cannot be claimed. Honest negative result.

**10.2/10.3/10.4 — what the pipeline actually implements (grep of plumber.R):** `RunEnrichr`/`EnrichrDotPlot` **0 calls** (enrichR attached, never used), `FindDMEs`/`FindMarkers` **0**, `p.adjust`/FDR **0**, `ModuleTraitCorrelation` **0**. PD-vs-control is **only** via `ModulePreservation`. ⇒ **W8**: the research question's downstream aims — pathway enrichment of the CACNA module (calcium/mito/lysosomal), quantitative PD-vs-HC differences — are **not answerable by the current pipeline**. Module preservation tests only whether control modules survive in PD; it yields no enrichment and no directional effect size. Multiple-testing (10.4) is moot — there is no per-module p-value stage.

**10.6 Negative controls — all candidate genes present + survive `fraction=0.05`:** SNCA 85.2%, PINK1 73.7%, PARK7 43.7%, GBA 30.8%, LRRK2 10.4%, PRKN 11.3%, and DA-identity TH 80.0% / SLC6A3 88.0% / SLC18A2 95.5%. ⇒ **N7**: no control is built in; recommend using the highly-expressed DA-identity trio (TH/SLC6A3/SLC18A2) as a positive control — they *should* co-cluster in a coherent module; if they don't, the network is untrustworthy. Their high detection (80–95%) also contrasts sharply with the target genes (CACNA1D 8%, CACNA1C 6%), reinforcing that the co-localization difficulty (B2) is driven by **target-gene sparsity** (W3), not a pipeline defect.

### Step 1 — §2 Environment & Dependency Audit

**W1 — Environment drift (PC/Docker vs verified env).** Native = R 4.5.3 + OpenBLAS + hdWGCNA 0.4.12 (tarball). Docker `r-service/Dockerfile` = `rocker/r-ver:4.4` + rocker default BLAS + `install_github('smorabit/hdWGCNA', ref='dev')` (unpinned HEAD). Different R minor, BLAS, and an unpinned hdWGCNA can shift module boundaries or behavior. This session's results (soft power 4; 7,465 genes) were on the native env and aren't guaranteed identical on the PC. *Fix (flag, ask first — reproducibility):* pin hdWGCNA to a dated commit, optionally pin R. Verify by fingerprinting the built image.

**W2 — Thread oversubscription.** `OMP_NUM_THREADS` unset ⇒ OpenBLAS threads run under `WGCNA_THREADS=4` × `future workers=2`. Perf/OOM risk on low-core machines. *User confident on PC hardware → ACCEPTED.* Cheap insurance if ever revisited: `OMP_NUM_THREADS=1` in compose, `WGCNA_THREADS` = PC core count.

**N1** — R 4.5.3 > tutorial/Docker 4.4; in-band (≥4.2), recorded for reproducibility.

### Step 2 — §3 Input Data Integrity Audit

Object: 28,437 genes × 8,216 cells; assay `RNA`; reductions `pca/harmony/umap`; layers `counts/data`. All cells labeled `Dopaminergic Neurons`.

**PIVOT CHECK — PASS.** `CACNA1D` (7.98%, 656 cells) and `CACNA1C` (6.18%, 508 cells) both present and confirmed in the 7,465-gene selected set. Project is answerable at defaults.

**W3 — `CACNA1C` thin margin.** Keep/drop by fraction (computed): kept ≤0.060, **dropped at 0.065**; falls out once fraction > ~0.062. `CACNA1D` drops at 0.10. Raising `fraction` even slightly (a plausible "be more stringent" move) silently removes Cav1.2 with no error. *Fix (flag — scientific param):* keep `fraction ≤ 0.06`; consider extending the existing 0-gene guard to assert target genes remain in `GetWGCNAGenes()`.

**W4 — PD underpowered.** 13 PD donors, only **4** clear `min_cells=100` (843/767/348/194); control 11 donors, **7** clear. Effective comparison 4 PD vs 7 control. Fragile to single-donor effects (pseudoreplication); biologically unavoidable (PD depletes DA neurons). Must be stated in the manuscript.

**N2** — `scale.data` dropped (CHANGELOG 0.3.0 memory fix). hdWGCNA scales metacells internally; verify no `ModuleEigengenes` breakage at §9.

### Step 3 — §4 Stage 1 `SetupForWGCNA`

**4.1 Parameters — PASS.** `gene_select`/`fraction`/`wgcna_name` wired from request with sensible defaults (`"fraction"`/0.05; plumber.R:128-129,241-245). Existing 0-gene tripwire at plumber.R:248-251. Target genes survive (see W3).

**4.2 Gene count — PASS (see N3).** 7,465 genes — inside the sane 1,000–10,000 band, on the higher end.

**4.3 No post-setup subsetting — PASS (clears BLOCKER risk).** Verified order in all 4 entrypoints: every `subset_condition()` runs **before** `SetupForWGCNA`.
- `run_full_pipeline`: setup only, no subset.
- `run_setup_through_soft_powers`: subset (489) → setup (492).
- `run_module_preservation`: subset ref+query (521-522) → setup ref (526) → `ProjectModules` query (534, internal setup). No `subset()` afterward.
- `run_gene_selection`: subset (424) → setup in loop on discarded `tmp` (438).

**N3** — 7,465 genes is on the high end (more noise / longer runtime than a tighter set). Judgment call, not a defect; the `/gene-selection` sweep exists to tune it. No action unless network quality (§8 grey-module fraction) looks poor.

### Step 4 — §5 Stage 2 `MetacellsByGroups` + `NormalizeMetacells`

Ran with the pipeline's exact params: `group.by=c("cell_type","sample_id")`, `reduction="pca"`, `k=25`, `max_shared=15`, `ident.group="cell_type"`, `min_cells=100` (default). Yield: **5,482 metacells** (PD 1,594 / control 3,888) from **11 donors** (4 PD, 7 control).

**5.1 group.by — PASS (clears #1 silent bug).** Donor is in `group.by`; all 5,482 metacells donor-pure (verified). No manufactured cross-patient co-expression.

**W5 — Metacell imbalance + donor concentration.** Per-donor metacell counts: PD `s.0098=684, s.0109=610, s.0116=229, s.0118=71`; control `s.0131=1000, s.0152=1000, s.0158=1000, s.0154=373, s.0159=210, s.0130=164, s.0142=141`. Three control donors hit the `target_metacells=1000` cap, so **3,000 of 3,888 control metacells come from just 3 donors**, and control outweighs PD 2.4:1. The pooled `/analyze` network — and the control reference network in `/module-preservation` — are dominated by a handful of high-cell-count donors. This compounds W4 (few donors) at the metacell level: co-expression estimates are effectively weighted toward those donors. *Not a blocker* (network is valid), but caps robustness and must temper any claim. Consider whether the primary result should come from `/module-preservation` (condition-specific) rather than the pooled `/analyze`.

**5.4 min_cells — surfaced, not silent.** 13 donor-groups dropped below 100 cells (hdWGCNA warned explicitly; also written to `donor_counts.csv`). Underpowering already tracked as W4.

**5.6 NormalizeMetacells — PASS.** Metacell `data` layer log-normalized (max 5.78); downstream `SetDatExpr(layer="data")` pulls it.

**N4 — Documented tutorial deviations (deliberate per CHANGELOG).** `reduction="pca"` not `harmony` — lower impact here because metacells are built *within* donor (KNN is within-donor, where cross-batch correction matters less). `max_shared=15` vs tutorial's 10 increases metacell overlap (more pseudo-replication among metacells). `target_metacells=1000` is hardcoded (not a request param), and is what drives the W5 caps. None are defects; flagged for awareness.

### Step 5 — §6 Stage 3 `SetDatExpr`

Ran the full chain (setup → metacells → normalize → `SetDatExpr`) and inspected the real `datExpr` via `GetDatExpr`.

**6.1 group selection — PASS.** `group_name="Dopaminergic Neurons"` (exact), `group.by="cell_type"` = the `ident.group` used in `MetacellsByGroups` (plumber.R:269-270). Tutorial's same-column requirement met.

**6.2 assay/layer — PASS.** `assay="RNA"`, `layer="data"` → normalized metacell matrix.

**6.3 multi-group intent — PASS.** Single `group_name` (not a vector) → network built on DA neurons only.

**6.4 matrix shape/quality — PASS.** `datExpr` = 5,482 metacells × 7,465 genes, correct WGCNA orientation (metacells in rows, genes in columns). Degenerate-gene scan: **0 all-NA, 0 any-NA, 0 all-zero, 0 zero-variance** — no gene will break `bicor`. Target genes present with real variance (`CACNA1D` var 0.040, in 4,286/5,482 metacells; `CACNA1C` var 0.035, in 3,896/5,482). Metacell aggregation lifted their detection from single-cell 8%/6% to 78%/71% — the correlation signal for the project's pivot genes is solid.

*No new findings — Section 6 is a clean pass.*

### Step 6 — §7 Stage 4 `TestSoftPowers` + soft-power selection

*Values verified earlier this session (fixed `set.seed(42)`, so reproducible); scratch artifacts since cleaned.*

**7.1 networkType consistency — PASS.** `p$network_type` (default `"signed"`) is the single source, used by both `TestSoftPowers` (plumber.R:318) and `ConstructNetwork` (plumber.R:361). No signed/unsigned split possible.

**7.2 power sweep ran — PASS.** `GetPowerTable` yields the full 20-row grid (`seq(1,10)` + `seq(12,30,2)`), written to `soft_power_table.csv`.

**7.3 chosen power justified — PASS.** Smallest power with SFT.R.sq ≥ 0.8 is **4** (R² = 0.837); `recommended_power=4`, `smallest_power=4`, `differ=false`. The recommendation reproduces `ConstructNetwork`'s own selection rule, so recommendation and built network agree. **Not hard-coded** — read from the table per dataset.

**7.4 reproducibility — PASS.** Power table + recommendation persisted to `out_dir` as `soft_power_table.csv`, `soft_power_recommendation.json`, `soft_power_plot.png`, `soft_power_plot.pdf`.

**N5 — soft power computed on the POOLED metacells.** `/test-soft-powers` without a condition subset builds the curve on all 5,482 metacells (PD+control). For `/module-preservation`, the reference network is built on ONE condition, whose scale-free curve can differ — so the power fed to `/module-preservation` should come from a `/test-soft-powers` call with `condition_col`+`ref_group` set (the service supports this). Using the pooled power for a condition-specific network is a subtle mismatch. Ties to W5 (the pooled network is control-dominated anyway).

**Verified-GOOD:** the no-power-reaches-0.8 case returns `null` + explicit warning, never `Inf` — directly resolves the audit's §7.3 concern that `ConstructNetwork(soft_power=NULL)` would otherwise silently auto-pick `Inf`.
