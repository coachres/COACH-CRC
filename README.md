# Integrative Multi-omics Analysis of Chinese Colorectal Cancers

## Introduction

This repository archives the computational workflows for the multi-omics
characterization of the COACH-CRC cohort (n = 1,005). The analysis code
encompasses pipelines for sequencing data processing, somatic mutation
detection, structural variation analysis, molecular subtyping, cis/trans
effect estimation, and prognostic modeling.

The analytical framework comprises the following components:

- **Primary data processing pipelines** for whole-genome sequencing (WGS) and
  whole-genome bisulfite sequencing (WGBS), together with somatic mutation and
  structural-variant calling.
- **Multi-omics consensus molecular subtyping (CMOS)** and an associated
  random forest classifier. CMOS is defined by NEMO consensus clustering of a
  five-omics matrix (transcriptome, focal copy number, DNA methylation,
  proteome, and phosphoproteome); the resulting classifier is trained on bulk
  RNA-seq data and applied to single-cell pseudobulk profiles.
- **Chromothripsis detection** based on structural variants and copy-number
  segments.
- **Cis and trans effect analysis** of somatic driver alterations on the
  multi-omics landscape, with statistical adjustment for tumor purity.
- **Survival and prognostic modeling**, including Cox regression and
  LASSO-based risk signatures with concordance-index (C-index) evaluation.
- **Publication figure generation**, including oncoprints (e.g., APC
  inversion), tumor-microenvironment (TME) heatmaps, and MSI panels
  integrating mutation, expression, and methylation data.

## Repository Organization

| Directory | Description |
|-----------|-------------|
| `molecular_subtyping/` | CMOS definition via MOVICS/NEMO and the associated random forest classifier. |
| `chromothripsis/` | ShatterSeek-based chromothripsis calling. |
| `cis_trans_analysis/` | Quantification of cis and trans effects of driver alterations. |
| `prognosis_analysis/` | Cox regression and prognostic model construction. |
| `plotting/` | Oncoprint, TME, and MSI/APC figure generation. |
| `pipeline/` | WGS / WGBS somatic and methylation processing; CNV/SV calling. |

## Module Descriptions

### Molecular subtyping and classification — `molecular_subtyping/`

- **`runMOVICS.R`** performs five-omics integration. Features are selected by
  median absolute deviation (retaining the top 10% for RNA, protein, and
  phosphosite features and the top 5% for methylation features), and the
  optimal cluster number is evaluated with `getClustNum`. A three-cluster NEMO
  solution yields the assigned subtypes **CMOS2–CMOS4**; samples lacking a
  confident assignment are labeled **CMOS1**. The script generates a
  multi-omics heatmap, a CMOS–CMS Sankey diagram, and Kaplan–Meier survival
  curves for overall survival (OS) and disease-free survival (DFS).
- **`CMOSclassifier.R`** trains a random forest model on bulk RNA-seq profiles
  (COACH-CRC) to discriminate the three assigned CMOS groups (CMOS2–CMOS4).
  Class imbalance is addressed through stratified subsampling, and predicted
  class probabilities are filtered using a confidence threshold to retain
  high-certainty predictions. The trained model is applied to single-cell
  pseudobulk profiles (microsatellite-stable samples) to assign CMOS labels
  for downstream immune-microenvironment comparison.

### Chromothripsis detection — `chromothripsis/`

- **`run_shatterseek.R`** wraps the **ShatterSeek** package (reference genome
  hg38). For each sample, the script ingests structural-variant and
  copy-number-segment tables, performs per-sample chromothripsis calling, and
  writes per-event statistics. Across samples, p-values for fragmentation
  joins, breakpoint enrichment, and chromosomal clustering are adjusted by the
  false discovery rate (FDR), and the results are merged into a single
  annotated table.

### Cis/trans effect analysis — `cis_trans_analysis/`

A linear-regression framework quantifies the associations between somatic
driver alterations and multi-omics features (transcriptome, proteome,
phosphoproteome, DNA methylation, focal CNV, and SV hotspots). Tumor purity is
included as a covariate in all models to distinguish somatic signal from
stromal contamination.

- **`01.snvindel.cis_effect.R`** — local (cis) effects of SNV/indel driver
  mutations.
- **`02.focal_cnv.cis_effect.R`** — cis effects of focal copy-number
  alterations.
- **`03.sv_hotspots.cis_effect.R`** — cis effects of structural-variant
  hotspots.
- **`04.trans_effect.R`** — distant (trans) effects of driver alterations
  across the same omics layers.

### Survival and prognostic modeling — `prognosis_analysis/`

- **`01.cox_regression.Rmd`** performs univariate and multivariate Cox
  proportional-hazards regression for OS and DFS, evaluated independently
  within each omics layer.
- **`02.prognostic_models_cindex.Rmd`** constructs and evaluates multi-omics
  prognostic signatures using LASSO regularization, reports concordance-index
  performance, and generates Kaplan–Meier survival plots.

### Visualization — `plotting/`

- **`plot_for_APC_inversion.R`** generates an APC oncoprint integrating
  coding, non-coding, SV, and CNV alterations with APC expression, compares
  APC expression across mutation classes, and renders a focused
  genomic-locus track of APC inversion hotspots.
- **`plotHeatmap_for_MSIsample.R`** assembles a multi-panel oncoprint and
  heatmap of mismatch-repair genes (mutation, protein, expression, and
  methylation) in microsatellite-instability (MSI) tumors, stratified by
  germline dMMR status, *MLH1* methylation, and *RAS*/*BRAF* status.
- **`plotHeatmap_for_CMOSfeature.R`** produces an integrated heatmap of CMOS
  subtypes with TME features, immune signatures (ESTIMATE / IOBR), pathway
  activities (GSVA Hallmark, PROGENy), and clinical and molecular annotations
  (CMS, WGD, TMB, SV count, epimutation burden, and driver mutations).

### Primary data processing — `pipeline/`

The Perl and Shell pipelines implement the primary processing of WGS and WGBS
data and downstream CNV/SV calling. They assume a conventional Linux layout
(described below) and are intended to be adapted to the local compute
environment.

- **`WGS_Pipeline_v1.0.pl`** — Sentieon-based WGS processing and somatic
  variant calling workflow.
- **`SomaticMutationDetect_v1.0.pl`** — somatic mutation detection and
  annotation workflow.
- **`WGBS_Pipeline_v1.0.pl`** — whole-genome bisulfite sequencing
  methylation-calling pipeline.
- **`runBRASS.sh`** — ascatNgs tumor/normal copy-number profiling followed by
  BRASS structural-variant detection.

## Software and Computational Environment

The workflows are designed for a standard Linux environment. External tools
and reference resources are expected under a single root: binaries under
`/opt/bioinfo` (with `/opt/bioinfo/database` for reference data), the R
runtime under `/opt/R`, and scratch space under `/data/scratch`. Path
variables declared at the top of each pipeline script should be edited to
match the target deployment.

### R (>= 4.1 recommended)

MOVICS, ShatterSeek, caret, randomForest, tidyverse, ComplexHeatmap,
survminer, survival, glmnet, ConsensusClusterPlus, circlize, igraph,
patchwork, ggsci, ggpubr, ggsurvplot, readxl, impute, extrafont, and their
dependencies, installed from CRAN and Bioconductor.

### External tools

- **WGS** (`WGS_Pipeline_v1.0.pl`, `SomaticMutationDetect_v1.0.pl`):
  Sentieon (requires a license server; set `SENTIEON_LICENSE`), SOAPnuke,
  GATK4, MuSE, Manta, Strelka2, Lancet, SvABA, SomaticSniper, bam-readcount,
  bedtools, vt, bcftools, samtools / bgzip / tabix, MSIsensor-pro, PCGR,
  BAM-matcher, and the CancerIT hts-ref bundle.
- **CNV / SV** (`runBRASS.sh`): ascatNgs, BRASS (CancerIT), and the CancerIT
  GRCh38 database bundle.
- **WGBS** (`WGBS_Pipeline_v1.0.pl`): BitMapperBS, MethylDackel, sambamba,
  Picard, samtools, and SOAPnuke.
- **Chromothripsis** (`chromothripsis/`): the **ShatterSeek** R package,
  together with structural-variant and copy-number-segment input tables.

The reference genome and GATK resource bundle (GRCh38 / hg38) and auxiliary
databases (COSMIC, PCGR data, microsatellite-site lists, region BED files,
etc.) are not distributed with this repository; these resources must be
obtained separately and referenced in the pipeline configuration.

## Data Preparation and Availability

Input data are referenced by file name within the analysis scripts and are
**not included** in this repository.

- `runMOVICS.R` expects its inputs (sample list, expression, focal CNV, DMR,
  protein, phosphosite, and clinical tables) under an `input/` subdirectory
  relative to the script.
- The scripts in `cis_trans_analysis/` and `prognosis_analysis/` reference
  intermediate and result tables by file name and expect them in the working
  directory (or via the path variables defined at the top of each script).
- Large input and output tables, `*.RData` intermediates, and session files
  are excluded from version control through `.gitignore`.

Processed data supporting the associated publication are available from the
authors upon reasonable request.

## Reproducibility Considerations

1. **Pipeline configuration.** The Perl/Shell pipelines in `pipeline/` are
   parameterized for tools under `/opt/bioinfo`, `/opt/R`, and
   `/data/scratch`; update the path variables at the top of each script when
   deploying to a different environment. `runBRASS.sh` uses
   `/path/to/...` placeholders for the two user-defined output directories.
2. **Input placement.** The analysis scripts resolve inputs by file name;
   ensure that the required tables are present in the working directory (or
   adjust the path variables) before execution. `runMOVICS.R` additionally
   requires an `input/` subdirectory.
3. **Cross-platform transfer.** The bulk-trained CMOS classifier is applied
   to single-cell pseudobulk profiles for *exploratory*
   immune-microenvironment comparison. Cross-platform batch effects and the
   absence of a single-cell gold standard should be taken into account when
   interpreting these assignments.

## Citation and License

If you use these workflows, please cite the associated COACH-CRC publication
(citation to be added). This repository is released under the MIT License;
the `LICENSE` file will be added prior to public release.
