# HRF Analysis of Resection-Labeled Clusters

This repository contains MATLAB and SLURM job scripts used to analyze fMRI activation clusters in epilepsy patients. Clusters are labeled based on their spatial relationship with the post-surgical resection mask, and hemodynamic response function (HRF) parameters are extracted from each labeled cluster.

## Overview

The analysis consists of two main stages:

### 1. Cluster Labeling

Each fMRI activation cluster is classified into one of the following categories:

- **Fully Concordant (FC)**: The voxel with the maximum Z-statistic lies *inside* the resection mask.
- **Partially Concordant (PC)**: The max voxel is *outside* the resection mask but within a 2 cm spherical dilation of it.
- **Partially Discordant (PD)**: The max voxel is outside, but some other part of the cluster overlaps with the resection mask.
- **Fully Discordant (FD)**: No part of the cluster intersects with the resection mask or its dilation.

### 2. HRF Feature Extraction

For each labeled cluster:

- HRFs are reconstructed voxel-wise using FLOBS basis functions (PE1, PE2, PE3).
- For each voxel:
  - Identify the **peak value** and **peak time**.
  - Compute the **FWHM** (full width at half maximum).
- Cluster-level summary includes:
  - **Mean HRF** across voxels
  - **Dominant peak** classification (positive or negative) using early vs late peak rules
  - **Standard deviation of peak latency**
  - **Standard deviation of FWHM**

## Repository Contents

- `*.m`: MATLAB scripts for HRF reconstruction and feature extraction.
- `*.slurm`: SLURM-compatible Bash scripts to run jobs on a computing cluster (e.g., ARC).

## Notes

- No fMRI or clinical data are included in this repository due to patient privacy.
- Subject IDs (e.g., `ICE070`) are anonymized and do not contain identifiable information.
- Paths are local to the user's ARC environment and should be adapted for reuse.



> **Note**: This project was conducted as part of a master’s research thesis in biomedical engineering, focusing on the use of fMRI for evaluating surgical outcomes in epilepsy.
