# Archived: human GSE243639 derived outputs

Archived 2026-08-13. Contains the **built/derived** artifacts from the human PD-vs-control
analysis, set aside so the pipeline can be rerun cleanly from raw data:

- `built_seurat_objects/` — `seurat_GSE243639_SNc.rds`, `seurat_GSE243639_SNc_clustered.rds`
  (moved from `DATA_UNZIPPED/`), plus `cluster_marker_scores.csv` and `markers_featureplot.png`
- `docker_shared_volume/` — `seurat_GSE243639_SNc.rds` (moved from `local_data/`, the Docker
  `/shared` mount)
- `analysis_results/` — the full prior `analysis_results/human_GSE243639/` tree (README,
  `pd_vs_control_combined_dme/`, `pd_vs_control_separate_networks/`)

**Not archived** — left in place in `DATA_UNZIPPED/` to rebuild from:
raw 10x counts (`SEURAT_INPUT/`, `GSE243639_RAW.tar`, `GSE243639_Filtered_count_table.csv.gz`),
clinical/metadata CSVs, and every `.R` build/analysis script
(`build_seurat_GSE243639.R`, `make_da_object.R`, `relabel_cluster7.R`,
`build_control_network.R`, `build_pd_network.R`, `run_dme_human_snc.R`).

Rerunning the human pipeline from scratch will regenerate everything in this folder at its
original paths; this archive is a snapshot of the pre-rebuild state, kept for comparison.
