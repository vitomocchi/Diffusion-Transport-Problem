# Diffusion-Transport Problem in MATLAB

## Introduction
This project tackles a diffusion-transport problem using MATLAB. It's divided into two parts: the first part focuses on calculating the exact solution and approximating it using three finite difference methods. The second part involves a variable substitution to compute exact approximated solutions.

## Part 1: Analysis
- **Exact Solution**: Obtained using MATLAB's symbolic toolbox.
- **Numerical Solution**: Approximated using finite difference methods.
- **Discretization Methods**: Centered, Forward, and Backward Differences.

## Code Summary
The main script `RunFD.m` sets up the problem parameters and computes numerical solutions using different finite difference schemes. It also plots the results and calculates the error norms.

## Results
Error norms for different schemes and parameters are tabulated, and graphical comparisons are provided.

## Conclusions
The project concludes with insights into the accuracy of different finite difference schemes, aligning with the concept of "upwind."
