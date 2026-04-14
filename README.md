# A Geometric Framework for 3D Leaf Movement by Orthonormal Bases: A Demonstration in Maranta leuconeura

## Scripts and Processed Data

This additional supporting information summarizes the processed datasets and analysis scripts used for the geometric analyses presented in the manuscript. These analyses operate **after** the computation of the orthonormal basis (ONB) from raw point-cloud data.

Raw point-cloud datasets are **not included** in this submission and **cannot be released**.  
High-resolution 3D point-cloud data may unintentionally reveal information about the experimental environment, acquisition setup, or other contextual details. Because the main conclusions of the study rely on the **geometric properties derived after ONB computation**, releasing the raw point clouds is neither required for reproducibility nor appropriate from a security and privacy standpoint.

All data necessary to reproduce the ONB-derived geometric analyses are fully provided below.

## Provided Materials

- **onb_analysis.py**  
  A custom Python module containing the utility functions used in the geometric analyses performed after ONB computation. This module is imported by the analysis notebooks.

- **20250713_df_vecs_only407.csv**  
  ONB-derived dataset used for the diurnal analysis.

- **20250905_all_df_summary.csv**  
  ONB-derived dataset used for the gravitropic-response analysis.

- **20251206_diurnal_reanalysis.ipynb**  
  Jupyter notebook for the diurnal analysis (input: 20250713_df_vecs_only407.csv).

- **20251205_reanalysis.ipynb**  
  Jupyter notebook for the gravitropic-response analysis (input: 20250905_all_df_summary.csv).

## Analysis Overview

- **Post-ONB geometric analysis**  
  Conducted using the processed datasets (20250713_df_vecs_only407.csv, 20250905_all_df_summary.csv) together with the helper functions implemented in onb_analysis.py.

- **Diurnal analysis**  
  - Input: 20250713_df_vecs_only407.csv  
  - Script: 20251206_diurnal_reanalysis.ipynb

- **Gravitropic response analysis**  
  - Input: 20250905_all_df_summary.csv  
  - Script: 20251205_reanalysis.ipynb

## Python Environment

The analyses were performed using the following environment:

Python 3.11.0  
matplotlib 3.8.4 (py311hca03da5_0)  
numpy 1.26.4 (py311he598dae_0)  
pandas 2.2.2 (py311h7aedaa7_0)  
scipy 1.13.1 (py311hac8794a_0)  
seaborn 0.13.2 (py311hca03da5_0)  
scikit-learn 1.4.2 (py311h7aedaa7_1)  
opencv 4.7.0 (py311ha1ab1f8_1)

## Notes on Reproducibility

The central claims of this study concern **geometric trajectories and orientation dynamics obtained after ONB computation**.  
These claims can be fully reproduced using the provided processed datasets and analysis scripts.  
The raw point-cloud data are not required to verify any of the scientific conclusions.

## Citation

A Geometric Framework for 3D Leaf Movement by Orthonormal Bases: A Demonstration in Maranta leuconeura

Miyuki T Nakata, Shotaro Sakita, Jion Shimoyama, Naoya Ando, Masahiro Takahara

Plant and Cell Physiology, pcag034 [https://doi.org/10.1093/pcp/pcag034](https://doi.org/10.1093/pcp/pcag034), 2026