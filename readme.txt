Title: Acoustic Hologram Optimisation Using Automatic Differentiation
Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
*: these authors contributed equally
Corresponding Author (email): Tatsuki Fushimi (tfushimi@slis.tsukuba.ac.jp)
Journal: Scientific Reports

Perquisite: 
Python (3.6.9): JAX, numpy, scipy, pandas
MATLAB (ver. 2020b): Image Processing Toolbox, k-wave (http://www.k-wave.org/), subtighplot(https://github.com/cortex-lab/MATLAB-tools/blob/master/burgbox/subtightplot.m) 
Julia (ver. 1.5.2): Adamopt, SpecialFunctions, Statistics, CSV, DataFrames, Random

Description:
This ZIP folder contains information necesarily to recreate the results in the manuscript:
- Diff-PAT Optimizer in Python (Automatic Differentiation)
- Diff-PAT Optimizer in Julia (Numerical Differentiation)
- Optimized performance data and execution time data for each optimizers
- Data analysis and figure plotting code for each figures in manuscript and supplementary material.

Some modifications were necessary to the code written by the authors of GS-PAT to extract the information (e.g. execution time) we needed, but those files are not included in this ZIP folder. This is due to the fact that we do not own the intellectual property over this code, and the licensing term for the codes are unclear. Readers wishing to obtain original codes for Eigensolver, Corrected Eigensolver, and GS-PAT are directed to the supplementary material of GS-PAT. Those who wish to see our C++ code for those optimizers are kindly requested to contact the corresponding author.  

The folders are put in the order in which it appears in the manuscript.
Step_1_PAL2AMP: Calculates the maximum theoretical pressure amplitude given a transducer array
Step_2_14x14x1 to step_4_32x32x1: Calculates the optimum field using Diff-PAT
Step_5_visualization: Program code to get necessarily information for discussion, and renders figures 2-5. 
Supplementary: Codes necessarily to obtain data for supplementary information are in these folders.

