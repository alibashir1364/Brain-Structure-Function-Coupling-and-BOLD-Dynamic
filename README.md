# Brain-Structure-Function-Coupling-and-BOLD-Dynamic
Structural Connectivity and Intrinsic BOLD Timescales

This repository contains the MATLAB code used to investigate how structural brain connectivity constrains the temporal dynamics of resting-state fMRI signals using a Graph Signal Processing (GSP) framework.

Overview

Structural connectivity provides the anatomical scaffold for large-scale brain dynamics. Previous studies have reported that brain regions with stronger structural connectivity tend to exhibit slower intrinsic functional dynamics. However, the signal components responsible for this relationship remain unclear.

In this project, resting-state BOLD activity is modeled as a graph signal supported on the structural connectome. Using Graph Fourier Transform (GFT)-based spectral decomposition, BOLD signals are separated into:

Structure-Coupled (Low-Frequency) Components
Structure-Decoupled (High-Frequency) Components

Intrinsic functional timescales are quantified using Relative Low-Frequency Power (RLFP) and related to regional structural connectivity strength.

Dataset

The analysis was performed using data from the Human Connectome Project (HCP):

100 unrelated healthy participants
Resting-state fMRI
Diffusion MRI (DWI)
Glasser 360-region cortical atlas
Main Analyses

The repository includes scripts for:

Construction of structural brain graphs
Graph Laplacian eigendecomposition
Graph Fourier Transform (GFT)
Spectral decomposition of BOLD signals
RLFP-based intrinsic timescale estimation
Structure–timescale correlation analysis
Subject-level statistical validation
Graph-based null model analysis
Transmodal vs. unimodal cortical comparisons
Publication-ready figure generation

Key Finding

The relationship between structural connectivity strength and intrinsic functional timescales is predominantly carried by structure-coupled (graph-smooth) BOLD activity. In contrast, structure-decoupled dynamics exhibit substantially weaker dependence on anatomical connectivity while retaining hierarchical organization across cortical systems.

Requirements
MATLAB R2022a or later
Statistics and Machine Learning Toolbox
Signal Processing Toolbox

Citation

If you use this code in your research, please cite:

Ali Bashirgonbadi, Mohamad Reza Salehi, and Hamid Soltanian-Zadeh, Senior Member, IEEE.
Structural Connectivity Selectively Constrains Intrinsic BOLD Timescales through Graph-Smooth Neural Activity.
(Manuscript under review)

License

This repository is released for academic and research purposes.
