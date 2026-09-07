# Prompt template — hdWGCNA Bridge work requests

Copy this template for any new request to Claude Code in this repo (new endpoint,
bugfix, refactor). Fill in each section before sending. Keep sections short —
the goal is precision, not length.

---

## Background

<!--
Context only — who/what/why. Reference CLAUDE.md instead of repeating it.
Example:
"This repo bridges hdWGCNA (R-only) to a Python backend via plumber, mirroring
the interface style of PyWGCNA where reasonable. I need to add support for
[X stage of the hdWGCNA pipeline]."
-->

## Request

<!--
One or two sentences, imperative verbs, scoped to a single change.
Example:
"Add a `/results/{job_id}/soft-power-plot` endpoint that returns the
TestSoftPowers() diagnostic plot as a PNG, following the same job-lookup
pattern as the existing /results/{job_id}/plot endpoint."
-->

## Inputs

<!--
List exact files, functions, or data this request depends on.
Example:
- r-service/plumber.R (existing /analyze endpoint, lines X-Y)
- hdWGCNA function: TestSoftPowers() — confirm signature before using it
- Sample metadata column name: "cell_type"
-->

## Deliverable

<!--
Format, length, language constraints.
Example:
"R code only, edits to r-service/plumber.R. No new files unless the existing
file would exceed ~150 lines. No explanation needed unless something is
ambiguous."
-->

## Guardrails

<!--
What NOT to do for this specific request, beyond what's already in CLAUDE.md.
Example:
"Do not change the /analyze endpoint's request/response shape. Do not add
caching. Do not add a new R package dependency without flagging it first."
-->

## Evaluation

<!--
Ask Claude to self-check before finalizing.
Example:
"Before finishing: confirm TestSoftPowers() returns a ggplot object compatible
with the same png()/print()/dev.off() pattern used elsewhere in plumber.R.
If you're not certain of its return type, say so instead of assuming."
-->
