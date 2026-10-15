# Multi-omics analysis pipeline for 1,005 Chinese colorectal cancers

Analysis code for multi-omics subtyping, molecular causality (cis/trans effects),
survival prognostic modeling, and figure generation in colorectal cancer (CRC),
developed for the COACH-CRC cohort.

## Introduction

This repository contains the R, Perl, and Shell scripts used to:

- Process raw **WGS and WGBS** data (Sentieon, somatic mutation detection, and
  methylation pipelines).
- Define a **five-omics molecular subtype (CMOS)** with
  [MOVICS](https://github.com/xlucpu/MOVICS) (NEMO consensus clustering), and
  classify samples—including single-cell pseudobulk profiles—with a random
  forest classifier trained on bulk transcriptomic data.
- Detect **chromothripsis** events (ShatterSeek) from structural variants and
  copy-number segments.
- Quantify **cis** and **trans** effects of somatic SNV/indel driver mutations
  , focal CNV, and SV hotspots on multi-omics features (RNA, protein,
  phosphosite, DNA methylation), controlling for tumor purity.
- Perform **survival analysis** (univariate and multivariate Cox regression for
  overall survival and disease-free survival) and construct prognostic models
  (LASSO + C-index) with Kaplan–Meier plots.

## Repository structure

| Directory | Contents |
|-----------|----------|
| `molecular_subtyping/` | `runMOVICS.R` — five-omics data integration and CMOS definition via MOVICS. `CMOSclassifier.R` — random forest classifier trained on bulk RNA, projected onto single-cell pseudobulk. 
| `pipeline/` | `WGS_Pipeline_v1.0.pl` (Sentieon-based WGS data process pipline), `SomaticMutationDetect_v1.0.pl`(somatic mutation detection pipeline), `WGBS_Pipeline_v1.0.pl` (WGBS methylation pipeline), `runBRASS.sh` (structural variant detection). |
| `shatterseek/` | `run_shatterseek.R` — chromothripsis detection with ShatterSeek across a sample list. |
| `cis_trans_analysis/` | `01.snvindel.cis_effect.R` (cis effects), `02.focal_cnv.cis_effect.R` (cis effects), `03.sv_hotspots.cis_effect.R` (cis effects); `04.trans_effect.R` (trans effects). Linear-regression framework with optional parallelization. |
| `prognosis_analysis/` | `01.cox_regression.Rmd` (univariate / multivariate Cox regression for OS and DFS); `02.prognostic_models_cindex.Rmd` (C-index evaluation, LASSO prognostic models, and Kaplan–Meier plots). |
| `plotting/` | `plot_for_APC_inversion.R` (APC inversion), `plotHeatmap_for_MSIsample.R` (mutation/expression/methylation landscape for MSI samples), `plotHeatmap_for_CMOSfeature.R` (CMOS tumor microenvironment / immune heatmap). |

### R (>= 4.1 recommended)

MOVICS, ShatterSeek, caret, randomForest, tidyverse, ComplexHeatmap, survminer,
survival, glmnet, ConsensusClusterPlus, circlize, igraph, patchwork, ggsci,
ggpubr, ggsurvplot, readxl, impute, extrafont, and related packages, installed
from CRAN and Bioconductor.

### Installation and dependencies

- **WGS** (`WGS_Pipeline_v1.0.pl`, `SomaticMutationDetect_v1.0.pl`): Sentieon
  (a license server is required; set `SENTIEON_LICENSE`), SOAPnuke, GATK4, MuSE,
  Manta, Strelka2, Lancet, SvABA, somatic-sniper, bam-readcount, bedtools, vt,
  bcftools, samtools, MSIsensor-pro, BAM-matcher, and the CancerIT hts-ref bundle.
- **SV** (`runBRASS.sh`): ascatNgs, BRASS (CancerIT), and the CancerIT GRCh38
  database bundle.
- **WGBS** (`WGBS_Pipeline_v1.0.pl`): BitMapperBS, MethylDackel, sambamba,
  picard, samtools, and SOAPnuke.
- **Chromothripsis** (`shatterseek/`): the **ShatterSeek** R package, together
  with structural-variant and copy-number-segment input files.

The reference genome / GATK resource bundle (GRCh38 / hg38) and the associated
database files (COSMIC, microsatellite sites, region BED files,etc.) are not 
distributed with this repository and must be obtained separately and referenced 
in the pipeline configuration.

## License

Distributed under the MIT License (see `LICENSE`).
