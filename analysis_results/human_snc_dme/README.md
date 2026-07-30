# Human SNc: Formal DME Test, PD vs Control (Combined Network)

Unlike the earlier `control_network/` and `pd_network/` folders (two *independently-built*
networks, not directly comparable module-for-module), this is **one combined network** built on
all 10,388 Dopaminergic Neurons pooled across both conditions (2,943 PD, 7,445 control), with
`FindDMEs` then comparing harmonized module eigengenes (hMEs) between the two conditions within
that single, fixed module definition — the statistically correct way to test "does this module's
activity differ by condition."

Efficiency note: earlier per-condition builds (`control_network/`, `pd_network/`) subset by
condition only, retaining all cell types (~80-90k cells each), and needed 25-34GB. Pre-filtering
to `cell_type="Dopaminergic Neurons"` first here cut this to ~10,388 cells total and only 6.07GB
peak — a large, avoidable cost in the earlier runs, now fixed for any future rebuilds.

## CACNA1D's module: turquoise, kME=0.53 — but NOT significantly differentially active

| Module | avg_log2FC | p.adj | Significant? |
|---|---|---|---|
| brown | −1.72 | 2.5×10⁻²² | Yes |
| blue | +0.87 | 5.0×10⁻⁴ | Yes |
| **turquoise (CACNA1D)** | −4.41 | **1.0** | **No** |

CACNA1D itself has strong connectivity within the turquoise module (kME=0.53, matching the
~0.54 seen in the earlier independently-built PD-only network), but the module's **overall
average activity does not differ significantly between PD and control** by this formal test.
Two *other* modules (brown, blue) do show significant differential activity.

This means the earlier "kME 0.147 (control) → 0.540 (PD)" finding — informative as evidence that
CACNA1D's co-expression neighborhood reorganizes across independently-built networks — is not
the same claim as "this module is differentially active," which the formal DME test does not
support here. See the top-level synthesis doc for how to talk about this distinction accurately.

Ignore the turquoise row's −4.41 log2FC as a real magnitude — module eigengenes can be negative
and near zero, making log-fold-change math unstable in that range (same caveat as the mouse
brown module in `mouse_GSE233866/snc_dme/`). The p-value, not this number, is what to trust.

## Files

`dendrogram.png`, `network_plot.png` (kME plot), `modules.csv` (combined network's gene→module
table), `dme_results.csv` (full DME table, all modules), `dme_volcano.png`, `dme_lollipop.png`,
soft-power diagnostics.
