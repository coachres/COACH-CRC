# Integrative Multi-omics Analysis of Chinese Colorectal Cancers

## Introduction

This repository archives the computational workflows for the multi-omics
characterization of the COACH-CRC cohort (n = 1,005). The analysis code
encompasses pipelines for sequencing data processing, somatic mutation
detection, structural variation analysis, molecular subtyping, cis/trans
effect estimation and prognostic modeling.

The analytical framework comprises the following components:

- **Data processing pipelines** for whole-genome sequencing (WGS) and
  whole-genome bisulfite sequencing (WGBS), together with somatic mutation and
  structural-variant calling.
- **Multi-omics consensus molecular subtyping (CMOS)** and an associated
  random forest classifier. CMOS is defined by NEMO consensus clustering of a
  five-omics matrix (transcriptome, focal copy number, DNA methylation,
  proteome, and phosphoproteome). The resulting classifier is trained on bulk
  RNA-seq data and applied to single-cell pseudobulk profiles.
- **Chromothripsis detection** based on structural variants and copy-number
  segments.
- **Cis and trans effect analysis** of somatic alterations on the multi-omics data.
- **Survival and prognostic modeling**, including Cox regression and
  LASSO-based risk signatures with concordance-index (C-index) evaluation.
- **Figure generation**, including oncoprints (e.g., APC inversion) and
  tumor-microenvironment (TME) heatmaps.

## Repository Organization

###  `molecular_subtyping/` : Molecular subtyping and classification

- **`runMOVICS.R`** performs five-omics integration. Features are selected by
  median absolute deviation (retaining the top 10% for RNA, protein, and
  phosphosite features and the top 5% for methylation features), and the
  optimal cluster number is evaluated with `getClustNum`. The script generates a
  multi-omics heatmap, a CMOS–CMS Sankey diagram, and Kaplan–Meier survival
  curves for overall survival (OS) and disease-free survival (DFS).
- **`CMOSclassifier.R`** trains a random forest model on bulk RNA-seq profiles
  (COACH-CRC) to discriminate the three assigned CMOS groups (CMOS2–CMOS4).
  The trained model is applied to single-cell pseudobulk profiles
  (for microsatellite-stable samples) to assign CMOS labels for downstream
  immune-microenvironment comparison.

###  `chromothripsis/` : Chromothripsis detection

- **`run_shatterseek.R`** performs chromothripsis calling with the ShatterSeek package
  (hg38) from structural variant and copy number segment tables.

###  `cis_trans_analysis/` : Cis/trans effect analysis

A linear-regression framework quantifies the associations between somatic
alterations and multi-omics features (transcriptome, proteome,
phosphoproteome, DNA methylation, focal CNV, and SV hotspots). 

- **`01.snvindel.cis_effect.R`** : cis effects of SNV/indel driver mutations.
- **`02.focal_cnv.cis_effect.R`** : cis effects of focal copy-number alterations.
- **`03.sv_hotspots.cis_effect.R`** : cis effects of structural-variant hotspots.
- **`04.trans_effect.R`** : trans effects of somatic alterations across the same omics data.

###  `prognosis_analysis/` : Survival and prognostic modeling

- **`01.cox_regression.Rmd`** performs univariate and multivariate Cox
  proportional-hazards regression for OS and DFS, evaluated independently
  within each omics layer.
- **`02.prognostic_models_cindex.Rmd`** constructs and evaluates multi-omics
  prognostic signatures using LASSO regularization, reports concordance-index
  performance, and generates Kaplan–Meier survival plots.

###  `plotting/` : Figure generation

- **`plot_for_APC_inversion.R`** : generates an APC oncoprint integrating
  coding, non-coding, SV, and CNV alterations with APC expression, compares
  APC expression across mutation classes, and renders a focused
  genomic-locus track of APC inversion hotspots.
- **`plotHeatmap_for_MSIsample.R`** : generates a multi-panel oncoprint and
  heatmap of mismatch-repair genes (mutation, protein, expression, and
  methylation) in microsatellite-instability (MSI) tumors, stratified by
  germline mismatch repair deficient (dMMR)  status, *MLH1* methylation,
  and *RAS*/*BRAF* status.
- **`plotHeatmap_for_CMOSfeature.R`** : generates an integrated heatmap of CMOS
  subtypes with TME features, immune signatures, pathway activities
  (GSVA Hallmark, PROGENy), and clinical and genomic features
  (CMS, WGD, TMB, SV count, epimutation burden, and driver mutations).

###  `pipeline/` : data processing and somatic variant calling

The Perl and Shell pipelines implement the primary processing of WGS and WGBS
data and downstream somatic variant calling.

- **`WGS_Pipeline_v1.0.pl`** : whole genome sequencing (WGS) data processing pipeline.
- **`SomaticMutationDetect_v1.0.pl`** : somatic variant calling workflow.
- **`WGBS_Pipeline_v1.0.pl`** : whole-genome bisulfite sequencing (WGBS) data processing and 
  methylation-calling pipeline.
- **`runBRASS.sh`** : structural variant detection by BRASS.

## Installation and dependencies

### R (>= 4.1 recommended)

MOVICS, ShatterSeek, caret, randomForest, tidyverse, ComplexHeatmap, survminer,
survival, glmnet, ConsensusClusterPlus, circlize, igraph, patchwork, ggsci,
ggpubr, ggsurvplot, readxl, impute, extrafont and related packages, installed
from CRAN and Bioconductor.

### External tools

The workflows are designed for a standard Linux environment and require the
following tools: Sentieon (a license server is required; set `SENTIEON_LICENSE`), 
SOAPnuke, GATK4, MuSE, Manta, Strelka2, Lancet, SvABA, somatic-sniper, 
bam-readcount, bedtools, bcftools, samtools, bgzip, tabix, MSIsensor-pro, 
BAM-matcher, ascatNgs, BRASS, BitMapperBS, MethylDackel, sambamba, and Picard.

The reference genome/GATK resource bundle (GRCh38/hg38) and associated
database files (COSMIC, microsatellite-site lists, region BED files, etc.) are
not distributed with this repository and must be obtained separately and
referenced in the pipeline configuration..

## Citation and License

If you use these workflows, please cite the associated COACH-CRC publication. 
This repository is released under the MIT License.
