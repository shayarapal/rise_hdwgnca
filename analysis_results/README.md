# CACNA1D Co-expression Network Analysis — Index

**Start here: [`SYNTHESIS.md`](SYNTHESIS.md)** — the cross-study synthesis explaining what each
analysis below shows, how strong its evidence is, and where it belongs in a paper.

## Folder layout

```
analysis_results/
├── SYNTHESIS.md                      ← read this first
├── human_GSE243639/                  Human SNc, PD vs. control (29 donors)
│   ├── pd_vs_control_separate_networks/    Two independently-built networks (qualitative)
│   └── pd_vs_control_combined_dme/         One combined network + formal DME test (rigorous)
├── mouse_GSE233866/                  Mouse midbrain, GSE233866
│   ├── healthy_baseline/             SNc vs. VTA, no disease
│   │   ├── snc_network/, vta_network/      Two independently-built networks
│   │   └── snc_vs_vta_preservation/        Formal preservation test
│   └── lesion_model/                 6-OHDA lesion (mouse PD model)
│       ├── snc_lesioned_network/, snc_intact_network/, vta_lesioned_network/
│       ├── lesioned_vs_intact_combined_dme_snc/     Disease vs. healthy, within SNc
│       └── snc_vs_vta_combined_dme_lesioned/        Region vs. region, during disease
└── reminders/
    └── reminder.md                   Running project log of decisions and pivots
```

Every subfolder with actual results has its own `README.md` with full detail (methods,
statistics, files, caveats). This index and `SYNTHESIS.md` are the only two documents meant to
be read top-to-bottom; everything else is reference material for a specific analysis.
