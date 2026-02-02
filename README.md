# TUD_Reasearch_Task_ATL
# Data-Driven Evaluation of Airport Slot Adherence

This repository contains the analysis code used to generate the results presented in the paper:

**“Data-Driven Evaluation of Airport Slot Adherence”**  
Tetiana Alifirenko, Rainer Koelle, Hartmut Fricke

The study proposes a data-driven framework for evaluating airport slot adherence using post-operations data, with a focus on empirically derived acceptable slot deviation windows and a hybrid buffer-based adherence indicator.

---

## Overview

Airport slot adherence analysis typically relies on fixed tolerance thresholds, which are not explicitly defined in post-operations datasets and vary across airports. This repository implements a reproducible analysis workflow to:

- quantify slot deviations using post-operations data,
- derive empirical and hybrid (empirical + Weibull-based) acceptable deviation windows,
- evaluate slot adherence under varying buffer definitions,
- compare adherence behaviour between CDM and non-CDM airports.

Due to data-sharing restrictions, raw operational data are **not included** in this repository.

---

## Repository Structure

```text
.
├── helpers.R              # Core helper functions (pre-processing, modeling, metrics)
├── analysis.R             # Main analysis pipeline
├── plots/                 # Generated figures used in the paper
├── tables/                # Generated tables used in the paper
├── sessionInfo.txt        # R session information for reproducibility
└── README.md
