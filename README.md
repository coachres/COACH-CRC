# Multi-omics analysis pipeline for 1,005 Chinese colorectal cancers

Analysis code for multi-omics subtyping, molecular causality (cis/trans effects),
survival prognostic modeling, and figure generation in colorectal cancer (CRC),
developed for the COACH-CRC cohort.

## Overview

This repository contains the R, Perl, and Shell scripts used to:

- Define a **five-omics molecular subtype (CMOS)** with
  [MOVICS](https://github.com/xlucpu/MOVICS) (NEMO consensus clustering), and
  classify samples—including single-cell pseudobulk profiles—with a random
  forest classifier trained on bulk transcriptomic data.
- Detect **chromothripsis** events (ShatterSeek) from structural variants and
  copy-number segments.
- Quantify **cis** and **trans** effects of somatic SNV/indel driver mutations
  on multi-omics features (RNA, protein, phosphosite, DNA methylation, focal
  CNV, and SV hotspots), controlling for tumor purity.
- Perform **survival analysis** (univariate and multivariate Cox regression for
  overall survival and disease-free survival) and construct prognostic models
  (LASSO + C-index) with Kaplan–Meier plots.
- Generate publication figures (oncoprints, tumor microenvironment / immune
  heatmaps, MSI panels).
- Process raw **WGS and WGBS** data (Sentieon, ascatNgs + BRASS, and
  methylation pipelines).

## Repository structure

| Directory | Contents |
|-----------|----------|
| `molecular_subtyping/` | `runMOVICS.R` — five-omics integration and CMOS definition via MOVICS/NEMO. `CMOSclassifier.R` — random forest classifier trained on bulk RNA, projected onto single-cell pseudobulk. `CMOSclassifier_v2.R` — revised version of the classifier with corrected preprocessing. |
| `pipeline/` | `WGS_Pipeline_v1.0.pl` (Sentieon-based WGS somatic pipeline), `SomaticMutationDetect_v1.0.pl`, `WGBS_Pipeline_v1.0.pl` (WGBS methylation pipeline), `runBRASS.sh` (ascatNgs CNV calling + BRASS structural variant detection). |
| `shatterseek/` | `run_shatterseek.R` — chromothripsis detection with ShatterSeek across a sample list. |
| `cis_trans_analysis/` | `01.snvindel.cis_effect.R`, `02.focal_cnv.cis_effect.R`, `03.sv_hotspots.cis_effect.R` (cis effects); `04.trans_effect.R` (trans effects). Linear-regression framework with optional parallelization. |
| `prognosis_analysis/` | `01.cox_regression.Rmd` (univariate / multivariate Cox regression for OS and DFS, per omics layer); `02.prognostic_models_cindex.Rmd` (C-index evaluation, LASSO prognostic models, Kaplan–Meier plots). |
| `plotting/` | `plot_for_APC_inversion.R` (APC inversion oncoprint), `plotHeatmap_for_MSIsample.R` (MSI mutation/expression/methylation panel), `plotHeatmap_for_CMOSfeature.R` (CMOS tumor microenvironment / immune heatmap). |

## Installation and dependencies

The scripts are designed for a standard Linux environment. Install the tools
listed below under a single root—by default `/opt/bioinfo` for binaries and
`/opt/bioinfo/database` for reference resources, with R under `/opt/R` and
scratch space under `/data/scratch`. The path variables at the top of each
pipeline script can be modified to match a different directory layout.

### R (>= 4.1 recommended)

MOVICS, ShatterSeek, caret, randomForest, tidyverse, ComplexHeatmap, survminer,
survival, glmnet, ConsensusClusterPlus, circlize, igraph, patchwork, ggsci,
ggpubr, ggsurvplot, readxl, impute, extrafont, and related packages, installed
from CRAN and Bioconductor.

### External tools

- **WGS** (`WGS_Pipeline_v1.0.pl`, `SomaticMutationDetect_v1.0.pl`): Sentieon
  (a license server is required; set `SENTIEON_LICENSE`), SOAPnuke, GATK4, MuSE,
  Manta, Strelka2, Lancet, SvABA, somatic-sniper, bam-readcount, bedtools, vt,
  bcftools, samtools / bgzip / tabix, MSIsensor-pro, PCGR, BAM-matcher, and the
  CancerIT hts-ref bundle.
- **CNV / SV** (`runBRASS.sh`): ascatNgs, BRASS (CancerIT), and the CancerIT
  GRCh38 database bundle.
- **WGBS** (`WGBS_Pipeline_v1.0.pl`): BitMapperBS, MethylDackel, sambamba,
  picard, samtools, and SOAPnuke.
- **Chromothripsis** (`shatterseek/`): the **ShatterSeek** R package, together
  with structural-variant and copy-number-segment input files.

The reference genome / GATK resource bundle (GRCh38 / hg38) and the associated
database files (COSMIC, PCGR data, microsatellite sites, region BED files,
etc.) are not distributed with this repository and must be obtained separately
and referenced in the pipeline configuration.

## Data availability

The scripts read external data files (`input/`, `00.input/`, `01.input/`,
expression, mutation, and clinical tables, and `*.RData` intermediates) that are
not included here. Processed data are available from the associated publication
and the corresponding authors upon reasonable request. Large input and output
files are excluded from the repository via `.gitignore`.

## Notes on reproducibility

- The Perl and Shell pipelines in `pipeline/` are configured for tools installed
  under `/opt/bioinfo`, `/opt/R`, and `/data/scratch`; update the path variables
  at the top of each script for a different environment. `runBRASS.sh` uses
  `/path/to/...` placeholders for the two user-defined output directories.
- `cis_trans_analysis/` and `prognosis_analysis/` require pre-built
  `00.input/*.RData` intermediates generated by the upstream processing steps.
- `CMOSclassifier_v2.R` supersedes `CMOSclassifier.R`; differences include
  consistent standardization of training, test, and single-cell data.

## License

Distributed under the MIT License (see `LICENSE`).
