# Analysis scripts

The per-study R scripts that produced the published results. They are **not** part of the
bridge service — the bridge is a general-purpose hdWGCNA wrapper, while these are the
specific analyses run against it for the CACNA1D SNc/VTA study.

They are tracked here so the code behind the paper is citable alongside the service.

## How they run

Each script begins:

```r
setwd("/app")
source("/app/plumber.R")
pipeline_libs()
```

`/app/plumber.R` is the bridge's own pipeline code, baked into the `r-service` image at
build time. The scripts call those functions directly in-process rather than over HTTP,
because a network build holds a multi-GB Seurat object and there is no reason to serialise
it across a socket.

Data is read by absolute container path (`/shared/DATA_GSE233866/...`,
`/shared/DATA_UNZIPPED/...`), which is the `SHARED_DIR` bind mount. Those two directory
names are load-bearing — renaming them breaks every script.

`docker-compose.yml` mounts only `SHARED_DIR`, so the container cannot see this directory.
To run a script, sync it into `SHARED_DIR` first:

```bash
make sync-analysis     # copies analysis/ -> $SHARED_DIR/DATA_*/
docker compose exec r-service Rscript /shared/DATA_GSE233866/build_snc_network.R
```

The copy in this repository is canonical. `make sync-analysis` is one-way — edit here, sync
out, never the reverse.

## mouse_GSE233866 — mouse midbrain snRNA-seq

Yaghmaeian Salmani et al. 6 animals (DatCre;TRAP, C57BL6, P90/P530), untreated baseline plus
a 6-OHDA lesioned/intact hemisphere arm.

| Order | Script | Purpose |
|---|---|---|
| 1 | `build_seurat_GSE233866.R` | Builds the untreated baseline object; splits SNc vs VTA |
| 1 | `build_seurat_lesion_intact.R` | Builds the lesioned-vs-intact hemisphere object |
| 2 | `diagnose_qc_floor.R` | Empirical re-derivation of the `nFeature_RNA` QC floor |
| 2 | `diagnose_cohort_counts.R` | Per-cohort cell counts |
| 2 | `diagnose_final_breakdown.R` | Final region/condition cell breakdown |
| 3 | `build_snc_network.R` | Healthy-baseline SNc co-expression network |
| 3 | `build_vta_network.R` | Healthy-baseline VTA network |
| 4 | `build_snc_intact_network.R` | Lesion-arm SNc network, intact hemisphere |
| 4 | `build_snc_lesioned_network.R` | Lesion-arm SNc network, lesioned hemisphere |
| 4 | `build_vta_intact_network.R` | Lesion-arm VTA network, intact hemisphere |
| 4 | `build_vta_lesioned_network.R` | Lesion-arm VTA network, lesioned hemisphere |
| 5 | `run_dme_snc.R` | DME, lesioned vs intact within SNc (one pooled network) |
| 5 | `run_dme_vta.R` | DME, lesioned vs intact within VTA |
| 5 | `run_dme_region_healthy.R` | DME, SNc vs VTA in the untreated cohort |
| 5 | `run_dme_region_intact.R` | DME, SNc vs VTA in the intact hemisphere |
| 5 | `run_dme_region_lesioned.R` | DME, SNc vs VTA in the lesioned hemisphere |
| 6 | `run_preservation_mouse.R` | Module preservation, SNc reference vs VTA query (200 permutations) |
| — | `extract_umap_combined.R` | Pulls per-cell UMAP coords, labels and CACNA1D expression for figures |
| — | `_scratch_inspect_tom.R` | Exploratory TOM inspection; not part of the published pipeline |

Scripts sharing an order number are independent of each other.

## human_GSE243639 — human SNc snRNA-seq

29 donors (15 Parkinson's, 14 control), 10x v3, 107,113 nuclei total.

| Order | Script | Purpose |
|---|---|---|
| 1 | `build_seurat_GSE243639.R` | Builds the full 29-donor object |
| 2 | `relabel_cluster7.R` | One-off relabel of cluster 7 to "Unassigned" |
| 3 | `make_da_object.R` | Subsets to the 8,216 dopaminergic nuclei the network is built on |
| 4 | `build_pd_network.R` | PD-donor network |
| 4 | `build_control_network.R` | Control-donor network |
| 5 | `run_dme_human_snc.R` | DME, PD vs control within one combined network |

## Outputs

Everything these scripts write — module tables, kME, DME statistics, preservation Z-scores,
enrichment exports and figures — lands in `$SHARED_DIR/analysis_results/`, outside this
repository. See the root README.
