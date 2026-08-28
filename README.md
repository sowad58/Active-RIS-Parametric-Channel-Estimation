# Parametric Channel Estimation and Design for Active-RIS-Assisted Communications

This repository contains the MATLAB simulation code for the paper:  
**"Parametric Channel Estimation and Design for Active-RIS-Assisted Communications"**  
*Authors: [Md. Shahriar Sadid], [A. A. Nasir], [S. Al-Ahmadi ],[S. Al-Ghadhban]*  
*Published in: IEEE Communications Letters, 2026*  
*Link to paper: [https://doi.org/10.1109/LCOMM.2026.3711925]*

## Overview
This code evaluates the performance of channel estimation and user tracking in an Active Reconfigurable Intelligent Surface (ARIS) empowered wireless network. Unlike passive RIS, Active RIS introduces dynamic thermal noise amplification, which severely degrades the performance of unstructured estimators (like Least Squares). 

This repository implements a rigorous **two-stage transmission protocol**:
1. **Offline Calibration (Stage 1):** A 2D-ESPRIT-based round-trip estimator to acquire the quasi-static BS-ARIS infrastructure prior without assuming perfect CSI.
2. **Online User Tracking (Stage 2):** A Parametric Maximum Likelihood Estimator (MLE) that directly extracts the physical geometric parameters (angles and distances) of the dominant Line-of-Sight (LoS) path in both near-field and far-field conditions.

The code generates the Achievable Rate and Normalized Mean Square Error (NMSE) comparisons against baseline unstructured estimators.

## File Structure
* **`main_code_v2.m`**: The main execution script. Run this file to perform the simulations and generate the Achievable Rate and NMSE plots.
* **`MLE3D_corrected.m`**: Implementation of the proposed near-field 3D Parametric MLE (estimates azimuth, elevation, and distance).
* **`MLE_corrected.m`**: Implementation of the far-field 2D Parametric MLE (estimates azimuth and elevation only).
* **`estimateBSARISChannel2DESPRIT.m`**: 2D-ESPRIT algorithm for the Stage 1 offline BS-ARIS channel calibration.
* **`generateBSARISCalibrationData.m`**: Generates the noisy calibration data for Stage 1.
* **`nearFieldChan.m` / `nearFieldChan_withAMplitude.m`**: Generates the exact near-field spherical wave channels based on the physical geometry.
* **`UPA_BasisElupnew.m` / `UPA_Codebook.m` / `UPA_Evaluation.m`**: Utility functions to generate the Uniform Planar Array (UPA) steering vectors and orthogonal codebooks.
* **`WideTwobeam32.mat`**: Pre-computed wide-beam initialization codebook to accelerate the spatial search.

## Requirements
* MATLAB (Tested on R2025)
* Parallel Computing Toolbox (optional, but recommended as the script utilizes `parfor` loops for faster Monte Carlo simulations).

## How to Run
1. Clone or download this repository to your local machine.
2. Ensure all `.m` files and the `.mat` file are in the same directory.
3. Open `main_code_v2.m` in MATLAB.
4. Set your desired simulation parameters at the top of the script (e.g., `NFConf`, `porpose`, `LSConf`).
5. Run the script. The script automatically averages and plots the results upon completion.

## Citation
If you use this code in your research, please cite our paper:
```bibtex
## Citation
If you use this code in your research, please cite our paper:

@ARTICLE{11602068,
  author={Sadid, Md. Shahriar and Nasir, Ali A. and Al-Ahmadi, Saad and Al-Ghadhban, Samir},
  journal={IEEE Communications Letters}, 
  title={Parametric Channel Estimation and Design for Active-RIS-Assisted Communications}, 
  year={2026},
  volume={30},
  pages={2595-2599},
  doi={10.1109/LCOMM.2026.3711925}
}
